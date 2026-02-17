% BRC06_OBJ3_PLSR_NestedCV_Tier1_Tier2_4PreprocCompare_v1.m
% -------------------------------------------------------------------------
% Objective 3 (public/audit-ready):
% Benchmark PLSR predictive performance under four spectral preprocessing
% pipelines using repeated nested cross-validation.
%
% Preprocessing configurations compared:
%   (1) RAW         : no preprocessing
%   (2) SNV         : Standard Normal Variate (row-wise)
%   (3) SNV_SG1st   : SNV + Savitzky–Golay 1st derivative (poly=2, window=11)
%   (4) SNV_SG2nd   : SNV + Savitzky–Golay 2nd derivative (poly=2, window=11)
%
% Validation design:
%   - Repeated nested CV:
%       Outer CV: K=5 folds, repeated R=50 times
%       Stratification: Maturity × N2 (to balance key experimental strata)
%       Optional guardrail: each outer TRAIN split must contain all Part levels
%   - Inner CV (LV selection): Kinner=4 folds; select LV by minimum RMSE
%     with LV upper bound: LVmax = min(15, nTrain-2, p)
%
% Outcomes (two tiers; edit lists below to match your matrix headers):
%   - Tier 1: primary endpoints (e.g., yield, phenolics, antioxidant assays)
%   - Tier 2: chemical class aggregates (SUM_* variables)
%
% Outputs (per tier and per preprocessing config):
%   - OBJ3_Summary_<Tier>_<Config>.xlsx
%       pooled + repeat-level summary metrics per response variable
%   - OBJ3_RepeatMetrics_<Tier>_<Config>.xlsx
%       metrics per repeat (R2, RMSE, MAE, Bias, RPD, LV_median)
%   - OBJ3_Predictions_<Tier>_<Config>.xlsx
%       long-format predictions; one Excel sheet per response variable Y
%   - OBJ3_VIPraw_<Tier>_<Config>.mat
%       raw VIP vectors per outer fit (for downstream VIP stability mapping)
%
% Additional tier-level comparison workbook:
%   - OBJ3_Compare_4CFG_<Tier>.xlsx
%       pooled/long summaries across all 4 preprocessing configurations
%
% Input requirements:
%   - Excel matrix with:
%       * spectral columns named nm_<integer> (e.g., nm_250 ... nm_1800)
%       * factor columns for Part, Maturity, N2 (names matched robustly)
%       * response columns for Tier 1 / Tier 2 variables
%
% MATLAB toolboxes:
%   - Statistics and Machine Learning Toolbox (plsregress, cvpartition)
%   - Signal Processing Toolbox (sgolay) for SG derivatives
%
% Notes:
%   - This script writes outputs to OUTDIR relative to the current working
%     directory (pwd). Run from your repo root or set OUTDIR explicitly.
%   - Default "core" filter (cultivar/extraction) can be edited below.
% -------------------------------------------------------------------------

clear; clc;

%% =========================
% USER PARAMETERS (EDIT)
% =========================
INPUT_XLSX  = 'Matriz_Brocoli_SUM_1nm_ASCII.xlsx';
SHEET_NAME  = 'Matriz';

OUTDIR      = fullfile(pwd, 'OBJ3_OUT_4CFG');
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

% Core definition (edit if needed)
CULTIVAR_CORE   = "Pathernon";
EXTRACTION_CORE = "Ultrasonido";

% Tier 1 (primary endpoints) — EXACT matrix headers
TIER1_LIST = string(["Extraction_yield","Total_phenolics","DPPH","ABTS","Antihypertensive_act."]);

% Tier 2 (class aggregates) — EXACT matrix headers
TIER2_LIST = string([
    "SUM_Amino_acids"
    "SUM_N_related"
    "SUM_OrgAcids"
    "SUM_GSL_aliphatic"
    "SUM_GSL_indolic"
    "SUM_GSL_breakdown"
    "SUM_Flavonols"
    "SUM_CQA"
    "SUM_Coumaroyl"
    "SUM_Sinapate_esters"
    "SUM_FA_saturated"
    "SUM_FA_unsaturated"
    "SUM_Oxylipins_oxygenated"
    "SUM_GSL_total"
    "SUM_Phenylpropanoids_total"
    "SUM_Lipids_total"
    "SUM_Phenolics_MS"
]);

% Y quality control
MIN_N_NONNAN = 20;
MAX_NAN_FRAC = 0.30;

% Outer CV (repeated)
K_OUTER   = 5;
R_REPEATS = 50;
BASE_SEED = 12345;

% Inner CV for LV selection
K_INNER    = 4;
LV_MAX_CAP = 15;
LV_POLICY  = "minRMSE";  %#ok<NASGU> (kept for audit; minRMSE implemented)

% Guardrail (recommended)
ENFORCE_PART_COVERAGE = true;
MAX_TRIES_FOLDS = 400;

% SG derivative parameters (SNV + SG derivatives)
SG_POLY_ORDER = 2;
SG_FRAME_LEN  = 11; % odd
DELTA_NM      = 1;  % 1 nm grid

% Standardise X after preprocessing (TRAIN stats)
SCALE_X = true;

%% =========================
% LOAD DATA (PRESERVE NAMES)
% =========================
fprintf('Loading: %s\n', INPUT_XLSX);
opts = detectImportOptions(INPUT_XLSX, 'Sheet', SHEET_NAME);
opts.VariableNamingRule = 'preserve';
T = readtable(INPUT_XLSX, opts);

fprintf('Rows: %d | Cols: %d\n', height(T), width(T));

%% =========================
% IDENTIFY FACTORS (ROBUST)
% =========================
colVariedad = pickVarName(T, ["Variedad","Variety","Cultivar"]);
colParte    = pickVarName(T, ["Parte","Part"]);
colMad      = pickVarName(T, ["Maduracion","Maduración","Maturity"]);
colN2       = pickVarName(T, ["Aplicacion_N2","Aplicación_N2","Aplicación N2","Aplicacion N2","N2","Nitrogen"]);
colExt      = pickVarName(T, ["Extraccion","Extracción","Extraction"]);

assert(colVariedad~="" && colExt~="" && colParte~="" && colMad~="" && colN2~="", ...
    'Missing required factor columns. Check matrix headers.');

%% =========================
% FILTER CORE
% =========================
varietyStr = string(T.(colVariedad));
extrStr    = string(T.(colExt));
isCore = strcmpi(strtrim(varietyStr), CULTIVAR_CORE) & strcmpi(strtrim(extrStr), EXTRACTION_CORE);
Tcore  = T(isCore, :);

fprintf('Core subset (Variedad=%s, Extraccion=%s): n = %d\n', CULTIVAR_CORE, EXTRACTION_CORE, height(Tcore));
assert(height(Tcore) >= 10, 'Core subset too small; check core filters.');

%% =========================
% EXTRACT X (SPECTRA) + WAVELENGTH VECTOR
% =========================
allNames  = string(Tcore.Properties.VariableNames);
isSpec    = startsWith(allNames, "nm_");
specNames = allNames(isSpec);
assert(~isempty(specNames), 'No spectral columns found (expected "nm_*").');

% Parse wavelengths and sort columns accordingly
wl = nan(numel(specNames),1);
for j=1:numel(specNames)
    tok = regexp(specNames(j), '^nm_(\d+)$', 'tokens', 'once');
    assert(~isempty(tok), 'Bad spectral name: %s', specNames(j));
    wl(j) = str2double(tok{1});
end
[wl, ord] = sort(wl,'ascend');
specNames = specNames(ord);

Xraw = double(table2array(Tcore(:, specNames)));

% Grid sanity
dx = median(diff(wl));
if any(abs(diff(wl) - dx) > 1e-9)
    warning('Wavelength grid not perfectly uniform. Using dx=%.6g for SG derivatives.', dx);
end
if abs(dx - DELTA_NM) > 1e-6
    warning('Detected dx=%.6g but DELTA_NM=%.6g. Consider setting DELTA_NM=dx.', dx, DELTA_NM);
end

%% =========================
% FACTORS & STRATA
% =========================
Part   = categorical(string(Tcore.(colParte)));
Mat    = categorical(string(Tcore.(colMad)));
N2     = categorical(string(Tcore.(colN2)));
strata = categorical(strcat(string(Mat), "×", string(N2)));

%% =========================
% CONFIGS (4-way comparison)
% =========================
configs = struct([]);

configs(1).name = "RAW";
configs(1).preprocessFcn = @(X) preprocess_raw(X);

configs(2).name = "SNV";
configs(2).preprocessFcn = @(X) preprocess_snv(X);

configs(3).name = "SNV_SG1st";
configs(3).preprocessFcn = @(X) preprocess_snv_sg(X, SG_POLY_ORDER, SG_FRAME_LEN, 1, DELTA_NM);

configs(4).name = "SNV_SG2nd";
configs(4).preprocessFcn = @(X) preprocess_snv_sg(X, SG_POLY_ORDER, SG_FRAME_LEN, 2, DELTA_NM);

%% =========================
% RUN BOTH TIERS
% =========================
runTier("Tier1", TIER1_LIST);
runTier("Tier2", TIER2_LIST);

fprintf('\nOBJ3 finished. Outputs in: %s\n', OUTDIR);

%% ============================================================
% LOCAL DRIVER
% ============================================================
function runTier(tierName, yList)
    % capture outer workspace vars
    OUTDIR   = evalin('base','OUTDIR');
    Tcore    = evalin('base','Tcore');
    Xraw     = evalin('base','Xraw');
    wl       = evalin('base','wl');
    Part     = evalin('base','Part');
    strata   = evalin('base','strata');
    configs  = evalin('base','configs');

    MIN_N_NONNAN = evalin('base','MIN_N_NONNAN');
    MAX_NAN_FRAC = evalin('base','MAX_NAN_FRAC');

    K_OUTER   = evalin('base','K_OUTER');
    R_REPEATS = evalin('base','R_REPEATS');
    BASE_SEED = evalin('base','BASE_SEED');
    K_INNER   = evalin('base','K_INNER');
    LV_MAX_CAP = evalin('base','LV_MAX_CAP');

    ENFORCE_PART_COVERAGE = evalin('base','ENFORCE_PART_COVERAGE');
    MAX_TRIES_FOLDS = evalin('base','MAX_TRIES_FOLDS');
    SCALE_X = evalin('base','SCALE_X');

    % Collect for tier-level comparison export
    allSummary = table();
    allRepeats = table();

    for c = 1:numel(configs)
        cfgName = string(configs(c).name);
        fprintf('\n=== OBJ3 | %s | %s ===\n', tierName, cfgName);

        outSummary = table('Size',[0 18], ...
    'VariableTypes', {'string','string','string','double', ...
                      'double','double','double','double','double', ...
                      'double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Tier','Config','Y','n', ...
                      'R2_pooled','RMSE_pooled','MAE_pooled','Bias_pooled','RPD_pooled', ...
                      'R2_mean','R2_sd','RMSE_mean','RMSE_sd','MAE_mean','MAE_sd','Bias_mean','Bias_sd','LV_median'});


        outRepeatsAll = table('Size',[0 10], ...
            'VariableTypes', {'string','string','string','double','double','double','double','double','double','double'}, ...
            'VariableNames', {'Tier','Config','Y','repeat','R2','RMSE','MAE','Bias','RPD','LV_median'});

        VIP_STORE = struct();

        predBook = fullfile(OUTDIR, sprintf('OBJ3_Predictions_%s_%s.xlsx', tierName, cfgName));
        if exist(predBook,'file'); delete(predBook); end

        for iY = 1:numel(yList)
            yName = yList(iY);

            if ~ismember(yName, string(Tcore.Properties.VariableNames))
                warning('%s Y not found in matrix: %s (skipping)', tierName, yName);
                continue;
            end

            yv = double(Tcore.(yName));
            nanFrac = mean(isnan(yv));
            nOK = sum(~isnan(yv));

            if nanFrac > MAX_NAN_FRAC || nOK < MIN_N_NONNAN || std(yv(~isnan(yv)))==0
                warning('Skipping Y=%s due to QC (nanFrac=%.2f, nOK=%d)', yName, nanFrac, nOK);
                continue;
            end

            valid = ~isnan(yv) & ~isundefined(Part) & ~isundefined(strata);
            X = Xraw(valid,:);
            y = yv(valid);
            part_v   = Part(valid);
            strata_v = strata(valid);

            params = struct();
            params.Kouter = K_OUTER;
            params.R      = R_REPEATS;
            params.baseSeed = BASE_SEED;
            params.Kinner = K_INNER;
            params.lvMaxCap = LV_MAX_CAP;
            params.enforcePartCoverage = ENFORCE_PART_COVERAGE;
            params.maxTriesFolds = MAX_TRIES_FOLDS;
            params.scaleX = SCALE_X;

            % Build ONE fold plan per Y (shared by all configs) => fair compare
            foldPlan = buildFoldPlan(strata_v, part_v, params);

            [predLong, repMetrics, pooledMetrics, vipAll, lvAll] = ...
                runNestedPLSR_givenFolds(X, y, foldPlan, configs(c).preprocessFcn, params);

            % Summary row
            newRow = { ...
                tierName, cfgName, yName, numel(y), ...
                pooledMetrics.R2, pooledMetrics.RMSE, pooledMetrics.MAE, pooledMetrics.Bias, pooledMetrics.RPD, ...
                mean(repMetrics.R2), std(repMetrics.R2), ...
                mean(repMetrics.RMSE), std(repMetrics.RMSE), ...
                mean(repMetrics.MAE), std(repMetrics.MAE), ...
                mean(repMetrics.Bias), std(repMetrics.Bias), ...
                median(lvAll) ...
                };

            outSummary = [outSummary; newRow]; %#ok<AGROW>

            % Repeats (add identifiers)
            repMetrics.Y = repmat(yName, height(repMetrics), 1);
            repMetrics.Tier = repmat(tierName, height(repMetrics), 1);
            repMetrics.Config = repmat(cfgName, height(repMetrics), 1);
            outRepeatsAll = [outRepeatsAll; repMetrics(:, {'Tier','Config','Y','repeat','R2','RMSE','MAE','Bias','RPD','LV_median'})]; %#ok<AGROW>

            % Predictions: one sheet per Y
            writetable(predLong, predBook, 'Sheet', safeSheet(yName));

            % VIP raw for downstream stability mapping
            yn = matlab.lang.makeValidName(yName);
            VIP_STORE.(yn).vip = vipAll;
            VIP_STORE.(yn).wl  = wl;
            VIP_STORE.(yn).lv  = lvAll;
            VIP_STORE.(yn).foldPlan = foldPlan; %#ok<STRNU> for audit
        end

        % Export per config
        writetable(outSummary, fullfile(OUTDIR, sprintf('OBJ3_Summary_%s_%s.xlsx', tierName, cfgName)), 'Sheet', 'Summary');
        writetable(outRepeatsAll, fullfile(OUTDIR, sprintf('OBJ3_RepeatMetrics_%s_%s.xlsx', tierName, cfgName)), 'Sheet', 'RepeatMetrics');
        save(fullfile(OUTDIR, sprintf('OBJ3_VIPraw_%s_%s.mat', tierName, cfgName)), 'VIP_STORE', '-v7.3');

        % Append to tier-level pools
        allSummary = [allSummary; outSummary]; %#ok<AGROW>
        allRepeats = [allRepeats; outRepeatsAll]; %#ok<AGROW>
    end

    % Tier-level comparison workbook
    cmpXlsx = fullfile(OUTDIR, sprintf('OBJ3_Compare_4CFG_%s.xlsx', tierName));
    if exist(cmpXlsx,'file'); delete(cmpXlsx); end
    writetable(allSummary, cmpXlsx, 'Sheet','Summary_AllCFG');
    writetable(allRepeats, cmpXlsx, 'Sheet','RepeatMetrics_AllCFG');

    % Optional wide-format quick view (pooled R2/RMSE by Y x Config)
    try
        TsubR2 = allSummary(:, {'Y','Config','R2_pooled'});
        TwR2 = unstack(TsubR2, 'R2_pooled', 'Config');
        writetable(TwR2, cmpXlsx, 'Sheet','R2pooled_wide');

        TsubRMSE = allSummary(:, {'Y','Config','RMSE_pooled'});
        TwRMSE = unstack(TsubRMSE, 'RMSE_pooled', 'Config');
        writetable(TwRMSE, cmpXlsx, 'Sheet','RMSEpooled_wide');
    catch
        warning('Wide-format export failed (unstack). Long format sheets are still available.');
    end

    fprintf('Tier compare workbook: %s\n', cmpXlsx);
end

function sh = safeSheet(name)
    name = char(name);
    name = regexprep(name, '[:\\/?*\[\]]', '_');
    if numel(name) > 31; sh = name(1:31); else; sh = name; end
end

%% ============================================================
% FOLD PLAN (identical across configs for each Y)
% ============================================================
function foldPlan = buildFoldPlan(strata, part, params)
    R = params.R;
    foldPlan = cell(R,1);

    for r = 1:R
        baseSeed = params.baseSeed + r;
        foldPlan{r} = makeStratifiedFolds_withCoverage(strata, part, params.Kouter, ...
            params.enforcePartCoverage, params.maxTriesFolds, baseSeed);
    end
end

function foldID = makeStratifiedFolds_withCoverage(strata, part, K, enforcePartCoverage, maxTries, seed0)
    n = numel(strata);
    foldID = zeros(n,1);
    partsAll = categories(part);

    for t = 1:maxTries
        rng(seed0 + 1000*t, 'twister');

        foldID(:) = 0;
        cats = categories(strata);

        for i=1:numel(cats)
            idx = find(strata == cats{i});
            idx = idx(randperm(numel(idx)));
            for j=1:numel(idx)
                foldID(idx(j)) = mod(j-1, K) + 1;
            end
        end

        if ~enforcePartCoverage
            return;
        end

        ok = true;
        for k=1:K
            trainIdx = (foldID ~= k);
            pTrain = categories(removecats(part(trainIdx)));
            if numel(pTrain) < numel(partsAll)
                ok = false; break;
            end
        end

        if ok
            return;
        end
    end

    error('Could not create stratified folds satisfying Part coverage after %d tries.', maxTries);
end

%% ============================================================
% NESTED PLSR (given foldPlan)
% ============================================================
function [predLong, repMetrics, pooledMetrics, vipAll, lvAll] = ...
    runNestedPLSR_givenFolds(X, y, foldPlan, preprocessFcn, params)

    n = numel(y);
    K = params.Kouter;
    R = params.R;

    predLong = table('Size',[0 6], ...
        'VariableTypes', {'double','double','double','double','double','double'}, ...
        'VariableNames', {'row','repeat','fold','y_true','y_pred','LV'});

    repMetrics = table('Size',[R 7], ...
        'VariableTypes', {'double','double','double','double','double','double','double'}, ...
        'VariableNames', {'repeat','R2','RMSE','MAE','Bias','RPD','LV_median'});

    vipAll = [];
    lvAll  = [];

    for r = 1:R
        foldID = foldPlan{r};
        yhat_r = nan(n,1);
        lv_r   = nan(n,1);

        for k = 1:K
            testIdx  = (foldID == k);
            trainIdx = ~testIdx;

            Xtr = X(trainIdx,:); ytr = y(trainIdx);
            Xte = X(testIdx,:);  yte = y(testIdx);

            % Preprocess (row-wise, no train stats)
            Xtr_p = preprocessFcn(Xtr);
            Xte_p = preprocessFcn(Xte);

            % Optional column scaling using TRAIN stats
            if params.scaleX
                [Xtr_p, muX, sdX] = centreScale(Xtr_p);
                Xte_p = applyCentreScale(Xte_p, muX, sdX);
            end

            % LV bounds
            lvMax = min(params.lvMaxCap, size(Xtr_p,1)-2);
            lvMax = max(1, min(lvMax, size(Xtr_p,2)));

            % Inner CV LV selection (min RMSE)
            lvSel = selectLV_innerCV(Xtr_p, ytr, params.Kinner, lvMax);

            % Fit + predict
            [~,~,~,~,beta,PCTVAR,~,stats] = plsregress(Xtr_p, ytr, lvSel);
            ypred = [ones(size(Xte_p,1),1), Xte_p] * beta;

            yhat_r(testIdx) = ypred;
            lv_r(testIdx)   = lvSel;

            % VIP (raw, per outer fit)
            vip = vip_plsr(stats, PCTVAR, size(Xtr_p,2));
            vipAll = [vipAll, vip(:)]; %#ok<AGROW>
            lvAll  = [lvAll; lvSel]; %#ok<AGROW>

            rowsTest = find(testIdx);
            tmp = table(rowsTest, repmat(r, numel(rowsTest),1), repmat(k, numel(rowsTest),1), yte, ypred, repmat(lvSel, numel(rowsTest),1), ...
                'VariableNames', {'row','repeat','fold','y_true','y_pred','LV'});
            predLong = [predLong; tmp]; %#ok<AGROW>
        end

        m = computeMetrics(y, yhat_r);
        repMetrics.repeat(r)    = r;
        repMetrics.R2(r)        = m.R2;
        repMetrics.RMSE(r)      = m.RMSE;
        repMetrics.MAE(r)       = m.MAE;
        repMetrics.Bias(r)      = m.Bias;
        repMetrics.RPD(r)       = m.RPD;
        repMetrics.LV_median(r) = median(lv_r, 'omitnan');
    end

    pooledMetrics = computeMetrics(predLong.y_true, predLong.y_pred);
end

%% ============================================================
% HELPERS
% ============================================================
function name = pickVarName(T, candidates)
    vars = string(T.Properties.VariableNames);
    name = "";
    for c = candidates
        idx = find(strcmpi(vars, c), 1);
        if ~isempty(idx); name = vars(idx); return; end
    end
    for c = candidates
        idx = find(contains(lower(vars), lower(c)), 1);
        if ~isempty(idx); name = vars(idx); return; end
    end
end

% --- Preprocessing
function Xout = preprocess_raw(X)
    Xout = X;
end

function Xsnv = preprocess_snv(X)
    mu = mean(X, 2, 'omitnan');
    sd = std(X, 0, 2, 'omitnan');
    sd(sd==0) = 1;
    Xsnv = (X - mu) ./ sd;
end

function Xout = preprocess_snv_sg(X, polyOrd, frameLen, derivOrd, delta)
    Xsnv = preprocess_snv(X);
    Xout = sg_derivative_rows_mirror(Xsnv, polyOrd, frameLen, derivOrd, delta);
end

function Xd = sg_derivative_rows_mirror(X, polyOrder, frameLen, derivOrder, dx)
    if mod(frameLen,2)==0
        error('SG_FRAME_LEN must be odd.');
    end
    if ~(derivOrder==1 || derivOrder==2)
        error('Only derivOrder=1 or 2 supported.');
    end

    [~, G] = sgolay(polyOrder, frameLen);
    h = factorial(derivOrder) / (dx^derivOrder) * G(:, derivOrder+1); % col vector

    half = (frameLen-1)/2;
    Xd = zeros(size(X));

    for i=1:size(X,1)
        x = X(i,:);
        left  = x(half:-1:1);
        right = x(end:-1:end-half+1);
        xpad = [left, x, right];

        y = conv(xpad, flipud(h), 'valid'); % length == numel(x)
        Xd(i,:) = y;
    end
end

% --- Scaling
function [Xcs, muX, sdX] = centreScale(X)
    muX = mean(X, 1, 'omitnan');
    sdX = std(X, 0, 1, 'omitnan');
    sdX(sdX==0) = 1;
    Xcs = (X - muX) ./ sdX;
end

function Xcs = applyCentreScale(X, muX, sdX)
    Xcs = (X - muX) ./ sdX;
end

% --- Inner CV LV selection (min RMSE)
function lvSel = selectLV_innerCV(Xtr, ytr, Kinner, lvMax)
    n = numel(ytr);
    cv = cvpartition(n, 'KFold', Kinner);

    rmse = nan(lvMax,1);
    for lv = 1:lvMax
        yhat = nan(n,1);
        for k=1:Kinner
            tr = training(cv,k);
            te = test(cv,k);
            [~,~,~,~,beta] = plsregress(Xtr(tr,:), ytr(tr), lv);
            yhat(te) = [ones(sum(te),1), Xtr(te,:)] * beta;
        end
        rmse(lv) = sqrt(mean((ytr - yhat).^2, 'omitnan'));
    end
    [~, lvSel] = min(rmse);
end

% --- VIP (raw; per outer fit)
function vip = vip_plsr(stats, PCTVAR, p)
    W = stats.W;               % p x LV
    LV = size(W,2);

    % Use Y variance explained per component (if available)
    if size(PCTVAR,1) >= 2
        yExp = PCTVAR(2,1:LV);
    else
        yExp = ones(1,LV)/LV;
    end

    denom = sum(yExp);
    if denom==0; denom=eps; end

    vip = zeros(p,1);
    for j=1:p
        vip(j) = sqrt( p * sum( yExp(:) .* (W(j,:).^2)' ) / denom );
    end
end

% --- Metrics
function m = computeMetrics(y, yhat)
    y = y(:); yhat = yhat(:);
    ok = ~isnan(y) & ~isnan(yhat);
    y = y(ok); yhat = yhat(ok);

    SSE = sum((y - yhat).^2);
    SST = sum((y - mean(y)).^2);
    m.R2   = 1 - SSE / max(SST, eps);
    m.RMSE = sqrt(mean((y - yhat).^2));
    m.MAE  = mean(abs(y - yhat));
    m.Bias = mean(yhat - y);
    m.RPD  = std(y) / max(m.RMSE, eps);
end
