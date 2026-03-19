function obj3_plsr_berry_composition()
% obj3_plsr_berry_composition
% ------------------------------------------------------------
% Objective 3 / Section 3.3
% Supervised PLS-R modelling of berry composition from
% hyperspectral predictors.
%
% This workflow:
%   1. reads the master Excel matrix and data dictionary,
%   2. extracts chemical endpoints,
%   3. builds predictor sets for FX10, FX17, and FUSION,
%   4. applies a preprocessing grid (RAW, SNV, SG1, SG2,
%      SNV_SG1, SNV_SG2),
%   5. evaluates models using repeated nested cross-validation,
%   6. applies outlier control within training folds only,
%   7. exports performance tables, predictions, figures, and a run log.
%
% Output location:
%   <pwd>/Objetivo_3/O3_PLSR_v22_<timestamp>/
%     - Polished_Performance.xlsx
%     - Polished_Predictions.xlsx
%     - Figures/
%     - run_log.txt
%
% Notes:
%   - White-background publication plots are enforced for PNG export.
%   - The analytical logic of the original workflow is preserved.
% ------------------------------------------------------------

%% ------------------------ CONFIG ---------------------------
cfg = struct();

% Input matrix
cfg.inFile    = fullfile(pwd, 'Matriz_CHEM_HSI_MASTER_96.xlsx');
cfg.sheetData = 'Matriz';
cfg.sheetDict = 'DataDictionary';

% Sensors and preprocessing grid
cfg.sensorsList = {'FX10', 'FX17', 'FUSION'};
cfg.preprocList = {'RAW', 'SNV', 'SG1', 'SG2', 'SNV_SG1', 'SNV_SG2'};

% Fusion cutoff (nm)
cfg.fusionCutoff = 950;

% Outlier settings
cfg.outlierSigma = 2.5;
cfg.maxDropPct   = 0.10;

% CV settings
cfg.cv.outerK = 10;
cfg.cv.outerR = 5;
cfg.cv.innerK = 10;
cfg.cv.seed   = 20260214;

% LV range
cfg.forceLV_Range = 1:12;

% Plot threshold
cfg.plotR2Threshold = 0.60;

% Output root (fixed folder name)
ts = datestr(now,'yyyymmdd_HHMMSS');
cfg.outBase = fullfile(pwd, 'Objetivo_3');
cfg.outDir  = fullfile(cfg.outBase, ['O3_PLSR_v22_' ts]);
cfg.figDir  = fullfile(cfg.outDir, 'Figures');

if ~exist(cfg.outBase,'dir'); mkdir(cfg.outBase); end
if ~exist(cfg.outDir,'dir');  mkdir(cfg.outDir);  end
if ~exist(cfg.figDir,'dir');  mkdir(cfg.figDir);  end

diary(fullfile(cfg.outDir,'run_log.txt'));

%% ---------------------- LOAD DATA --------------------------
fprintf('--- O3 v22 FINAL POLISHED RUN START ---\n');
fprintf('Input:  %s\n', cfg.inFile);
fprintf('Output: %s\n', cfg.outDir);

if ~isfile(cfg.inFile)
    diary off;
    error('Input file not found: %s', cfg.inFile);
end

opts = detectImportOptions(cfg.inFile, 'Sheet', cfg.sheetData);
opts.VariableNamingRule = 'preserve';
T = readtable(cfg.inFile, opts, 'Sheet', cfg.sheetData);

optsDD = detectImportOptions(cfg.inFile, 'Sheet', cfg.sheetDict);
optsDD.VariableNamingRule = 'preserve';
DD = readtable(cfg.inFile, optsDD, 'Sheet', cfg.sheetDict);

% Endpoints: Category == "Chem endpoint"
isChem = strcmpi(string(DD.Category), "Chem endpoint");
endpointsAll = string(DD.Variable(isChem));

vars = string(T.Properties.VariableNames);

%% ------------------- PREPARE SPECTRA -----------------------
% FX10
cols10 = startsWith(vars, "FX10_nm_");
wl10 = parseWavelengths(vars(cols10), "FX10");
[wl10, idx10] = sort(wl10);
X10_raw_all = table2array(T(:, cols10));
X10_raw_all = X10_raw_all(:, idx10);
dt10 = median(diff(wl10));

% FX17
cols17 = startsWith(vars, "FX17_nm_");
wl17 = parseWavelengths(vars(cols17), "FX17");
[wl17, idx17] = sort(wl17);
X17_raw_all = table2array(T(:, cols17));
X17_raw_all = X17_raw_all(:, idx17);
dt17 = median(diff(wl17));

% FUSION (concatenate RAW)
keep10 = wl10 <= cfg.fusionCutoff;
keep17 = wl17 > cfg.fusionCutoff;

X10_cut = X10_raw_all(:, keep10);
X17_cut = X17_raw_all(:, keep17);

X_fused_raw_all = [X10_cut, X17_cut];
n10_cut = size(X10_cut, 2);

dt10_cut = median(diff(wl10(keep10)));
dt17_cut = median(diff(wl17(keep17)));

%% -------------------- MAIN LOOP -----------------------
PerfRows = {};
PredRows = {};

for e = 1:numel(endpointsAll)
    endpoint = endpointsAll(e);
    if ~any(strcmp(vars, endpoint))
        continue;
    end

    y_full = T.(endpoint);
    valid = ~isnan(y_full);
    if sum(valid) < 20
        continue;
    end

    y = y_full(valid);
    fprintf('\n>>> MODELING: %s (n=%d) <<<\n', endpoint, numel(y));

    X10_curr  = X10_raw_all(valid, :);
    X17_curr  = X17_raw_all(valid, :);
    XFus_curr = X_fused_raw_all(valid, :);
    origIndices = find(valid);

    for s = 1:numel(cfg.sensorsList)
        sensorName = cfg.sensorsList{s};

        for p = 1:numel(cfg.preprocList)
            preproc = cfg.preprocList{p};

            % 1) PREPROCESS
            try
                if strcmp(sensorName, 'FX10')
                    X_pre = applyPolishedPreproc(X10_curr, preproc, dt10, false);
                elseif strcmp(sensorName, 'FX17')
                    X_pre = applyPolishedPreproc(X17_curr, preproc, dt17, false);
                else % FUSION
                    X_pre = applyFusionPreproc(XFus_curr, preproc, n10_cut, dt10_cut, dt17_cut);
                end
            catch ME
                fprintf('Preproc Error [%s|%s]: %s\n', sensorName, preproc, ME.message);
                continue;
            end

            % 2) CV
            cvOut = runPolishedCV(X_pre, y, cfg);

            % 3) STORE PERFORMANCE
            perfRow = {char(endpoint), sensorName, preproc, numel(y), ...
                       cvOut.R2cv, cvOut.RMSEcv, cvOut.RPD, cvOut.Bias, ...
                       cvOut.LV_Mean, cvOut.LV_Median, cvOut.AvgOutliers};
            PerfRows(end+1, :) = perfRow; %#ok<AGROW>

            % 4) STORE PREDICTIONS
            for i = 1:numel(y)
                PredRows{end+1, 1} = {char(endpoint), sensorName, preproc, ...
                    origIndices(i), y(i), cvOut.yhatMean(i)}; %#ok<AGROW>
            end

            % 5) PLOT (white-safe PNG)
            if cvOut.R2cv > cfg.plotR2Threshold
                fprintf('  + %s | %s: R2cv=%.2f (AvgOut=%.1f)\n', sensorName, preproc, cvOut.R2cv, cvOut.AvgOutliers);
                makePlot_Q1(y, cvOut.yhatMean, char(endpoint), sensorName, preproc, cvOut.R2cv, cfg.figDir);
            end
        end
    end
end

%% -------------------- EXPORT -----------------------
if isempty(PerfRows)
    diary off;
    error('No models succeeded. Check input matrix and dictionary.');
end

PerfT = cell2table(PerfRows, 'VariableNames', ...
    {'Endpoint','Sensor','Preprocess','nSamples','R2cv','RMSEcv','RPD','Bias','LV_Mean','LV_Median','AvgOutliersRemoved'});
writetable(PerfT, fullfile(cfg.outDir, 'Polished_Performance.xlsx'));

PredsT = vertcat(PredRows{:});
PredsT = cell2table(PredsT, 'VariableNames', ...
    {'Endpoint','Sensor','Preprocess','OriginalRowIdx','Reference','Prediction'});
writetable(PredsT, fullfile(cfg.outDir, 'Polished_Predictions.xlsx'));

fprintf('\nDone. Results saved in:\n  %s\n', cfg.outDir);
diary off;

end

%% ========================= LOCALS ===========================

function out = runPolishedCV(X, y, cfg)
n = numel(y);
yhat = nan(n, cfg.cv.outerR);
lvSelected = nan(cfg.cv.outerR, cfg.cv.outerK);
nOutliersTotal = 0;

for r = 1:cfg.cv.outerR
    rng(cfg.cv.seed + r);
    cvOuter = cvpartition(n, 'KFold', cfg.cv.outerK);

    for k = 1:cfg.cv.outerK
        idxTr = training(cvOuter, k);
        idxTe = test(cvOuter, k);

        Xtr = X(idxTr, :); ytr = y(idxTr);
        Xte = X(idxTe, :);

        % Outlier detection
        currentSeed = cfg.cv.seed + 10000*r + k;
        [cleanIdx, ~] = detectOutliersCV(Xtr, ytr, cfg.outlierSigma, cfg.maxDropPct, currentSeed);

        Xtr_clean = Xtr(cleanIdx, :);
        ytr_clean = ytr(cleanIdx);
        nOutliersTotal = nOutliersTotal + sum(~cleanIdx);

        % LV selection (strict one-SE fairness)
        bestLV = selectLV_OneSE_Strict(Xtr_clean, ytr_clean, cfg.cv.innerK, cfg.forceLV_Range, currentSeed);

        % Train & predict
        [~,~,~,~,beta] = plsregress(Xtr_clean, ytr_clean, bestLV);
        yhat(idxTe, r) = [ones(sum(idxTe),1) Xte] * beta;
        lvSelected(r, k) = bestLV;
    end
end

yhatMean = mean(yhat, 2, 'omitnan');
out.yhatMean = yhatMean;
out.RMSEcv = sqrt(mean((yhatMean - y).^2, 'omitnan'));
out.R2cv   = 1 - sum((yhatMean - y).^2, 'omitnan') / sum((y - mean(y,'omitnan')).^2, 'omitnan');
out.RPD    = std(y, 'omitnan')/out.RMSEcv;
out.Bias   = mean(yhatMean - y, 'omitnan');
out.LV_Mean   = mean(lvSelected(:), 'omitnan');
out.LV_Median = median(lvSelected(:), 'omitnan');
out.AvgOutliers = nOutliersTotal / (cfg.cv.outerR * cfg.cv.outerK);
end

function [cleanIdx, isOut] = detectOutliersCV(X, y, sigmaThresh, maxDropPct, seed)
if nargin >= 5; rng(seed); end

n = numel(y);
if n < 10
    cleanIdx = true(n,1); isOut = false(n,1); return;
end

K = 5;
cv = cvpartition(n, 'KFold', K);
yhat = nan(n,1);

globalMaxPilot = 6;

for k = 1:K
    tr = training(cv, k); te = test(cv, k);

    pilotLV_fold = min([globalMaxPilot, sum(tr)-1, rank(X(tr,:))]);
    pilotLV_fold = max(pilotLV_fold, 1);

    [~,~,~,~,beta] = plsregress(X(tr,:), y(tr), pilotLV_fold);
    yhat(te) = [ones(sum(te),1) X(te,:)] * beta;
end

resid = y - yhat;

medResid = median(resid, 'omitnan');
madVal = median(abs(resid - medResid), 'omitnan');
robustSigma = 1.4826 * madVal;

if ~isfinite(robustSigma) || robustSigma < eps
    robustSigma = std(resid, 'omitnan');
    if ~isfinite(robustSigma) || robustSigma < eps
        robustSigma = 1;
    end
end

zScore = abs(resid - medResid) / robustSigma;
isOut = zScore > sigmaThresh;

if sum(isOut) > floor(n * maxDropPct)
    [~, sortIdx] = sort(zScore, 'descend');
    maxAllowed = floor(n * maxDropPct);
    isOut(:) = false;
    isOut(sortIdx(1:maxAllowed)) = true;
end

cleanIdx = ~isOut;
end

function bestLV = selectLV_OneSE_Strict(X, y, kFold, lvRange, seed)
n = size(X,1);
rng(seed);
cv = cvpartition(n, 'KFold', kFold);

lvRange = lvRange(lvRange >= 1);
if isempty(lvRange); bestLV = 1; return; end

mse_folds = nan(cv.NumTestSets, max(lvRange));

for k = 1:cv.NumTestSets
    tr = training(cv, k); te = test(cv, k);

    maxLV_Fold = min([max(lvRange), sum(tr)-1, size(X,2), rank(X(tr,:))]);
    lvEff = lvRange(lvRange <= maxLV_Fold);
    if isempty(lvEff); lvEff = 1; end

    for iLV = 1:numel(lvEff)
        lv = lvEff(iLV);
        [~,~,~,~,beta] = plsregress(X(tr,:), y(tr), lv);
        yhat = [ones(sum(te),1) X(te,:)] * beta;
        mse_folds(k, lv) = mean((y(te) - yhat).^2, 'omitnan');
    end
end

validCols = find(all(~isnan(mse_folds), 1));
if isempty(validCols)
    bestLV = 1; return;
end

mse_folds = mse_folds(:, validCols);
lv_axis   = validCols;

mse_mean = mean(mse_folds, 1, 'omitnan');
mse_se   = std(mse_folds, 0, 1, 'omitnan') / sqrt(cv.NumTestSets);

[minMSE, idxMin] = min(mse_mean);
target = minMSE + mse_se(idxMin);

idxBest = find(mse_mean <= target, 1, 'first');
bestLV = lv_axis(idxBest);
end

function Xp = applyPolishedPreproc(X, preproc, dt, isFusionBlock)
window = 15;
order  = 2;
if nargin < 3 || isempty(dt); dt = 1; end

switch upper(preproc)
    case 'RAW'
        Xp = X;
    case 'SNV'
        Xp = doSNV(X);
    case 'SG1'
        Xp = doSG_Deriv_Padded(X, order, window, 1, dt);
    case 'SG2'
        Xp = doSG_Deriv_Padded(X, order, window, 2, dt);
    case 'SNV_SG1'
        Xtemp = doSNV(X);
        Xp = doSG_Deriv_Padded(Xtemp, order, window, 1, dt);
    case 'SNV_SG2'
        Xtemp = doSNV(X);
        Xp = doSG_Deriv_Padded(Xtemp, order, window, 2, dt);
    otherwise
        error('Unknown preproc: %s', preproc);
end

if isFusionBlock
    Xp = rowL2norm(Xp);
end
end

function Xp = applyFusionPreproc(Xfus, preproc, n10, dt10, dt17)
X10 = Xfus(:, 1:n10);
X17 = Xfus(:, n10+1:end);

switch upper(preproc)
    case 'RAW'
        Xp = [X10, X17];
        Xp = rowL2norm(Xp);

    case 'SNV'
        X10p = doSNV(X10); X10p = rowL2norm(X10p);
        X17p = doSNV(X17); X17p = rowL2norm(X17p);
        Xp = [X10p, X17p];

    case {'SG1','SG2','SNV_SG1','SNV_SG2'}
        X10p = applyPolishedPreproc(X10, preproc, dt10, true);
        X17p = applyPolishedPreproc(X17, preproc, dt17, true);
        Xp = [X10p, X17p];

    otherwise
        error('Unknown fusion preproc: %s', preproc);
end
end

function Xd = doSG_Deriv_Padded(X, order, window, deriv, dt)
L = size(X,2);
if L < window
    error('SG derivative: spectrum length (%d) < window (%d).', L, window);
end

[~, g] = sgolay(order, window);
scale  = factorial(deriv) / ((-dt)^deriv);
kernel = flipud(g(:, deriv + 1) * scale);

half = floor(window/2);
Xd = zeros(size(X));

for i = 1:size(X,1)
    row = X(i,:);
    left  = row(half+1:-1:2);
    right = row(end-1:-1:end-half);
    rowPad = [left, row, right];
    Xd(i,:) = conv(rowPad, kernel, 'valid');
end
end

function Xsnv = doSNV(X)
mu = mean(X, 2, 'omitnan');
sd = std(X, 0, 2, 'omitnan');
sd(sd < 1e-12) = 1;
Xsnv = (X - mu) ./ sd;
end

function Xb = rowL2norm(Xb)
nf = vecnorm(Xb, 2, 2);
nf(nf < eps) = 1;
Xb = Xb ./ nf;
end

function wl = parseWavelengths(varNames, sensor)
prefix = sensor + "_nm_";
wl = nan(size(varNames));
for i = 1:numel(varNames)
    s = erase(varNames(i), prefix);
    s = replace(s, "p", ".");
    s = replace(s, "_", ".");
    wl(i) = str2double(s);
end
end

function makePlot_Q1(y, yhat, ep, sensor, preproc, r2, outDir)
% Q1-style Pred vs Ref plot (WHITE background, dark-safe) -> .png + .fig
% Forces light theme regardless of MATLAB UI theme.

if ~exist(outDir,'dir'); mkdir(outDir); end

% Figure + axes forced to white
f = figure('Visible','off','Color','w','Units','pixels','Position',[100 100 650 520]);
set(f,'InvertHardcopy','off');

ax = axes(f); %#ok<LAXES>
set(ax,'Color','w','Box','on','LineWidth',1.1,'TickDir','out');
set(ax,'FontName','Times New Roman','FontSize',13,'FontWeight','bold');
set(ax,'XColor','k','YColor','k');
set(ax,'XGrid','on','YGrid','on','GridAlpha',1,'GridColor',[0.85 0.85 0.85]);

scatter(ax, y, yhat, 42, 'filled', ...
    'MarkerFaceColor',[0.15 0.55 0.25], 'MarkerEdgeColor','k', 'LineWidth',0.7);
hold(ax,'on');

mn = min([y(:); yhat(:)]);
mx = max([y(:); yhat(:)]);
plot(ax, [mn mx], [mn mx], 'k--', 'LineWidth', 1.2);

xlabel(ax,'Reference','FontName','Times New Roman','FontSize',16,'FontWeight','bold','Color','k');
ylabel(ax,'Prediction','FontName','Times New Roman','FontSize',16,'FontWeight','bold','Color','k');

title(ax, sprintf('%s | %s | %s (R2cv=%.2f)', ep, sensor, preproc, r2), ...
    'Interpreter','none','FontName','Times New Roman','FontSize',14,'FontWeight','bold','Color','k');

pad = 0.03*(mx-mn);
xlim(ax, [mn-pad, mx+pad]);
ylim(ax, [mn-pad, mx+pad]);

fnameBase = sprintf('Plot_%s_%s_%s', ep, sensor, preproc);
fnameBase = regexprep(fnameBase,'[^a-zA-Z0-9_\-]','_');

figPath = fullfile(outDir, [fnameBase '.fig']);
pngPath = fullfile(outDir, [fnameBase '.png']);

savefig(f, figPath);

try
    exportgraphics(f, pngPath, 'Resolution', 400, 'BackgroundColor','white');
catch
    set(f,'PaperPositionMode','auto');
    print(f, pngPath, '-dpng', '-r400');
end

close(f);
end