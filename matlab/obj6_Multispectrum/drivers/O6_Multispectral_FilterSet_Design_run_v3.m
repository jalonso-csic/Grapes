function O6_Multispectral_FilterSet_Design_run_v3()
% O6_Multispectral_FilterSet_Design_run_v3
% -------------------------------------------------------------------------
% Objective 6: Translate VIP interpretability into deployable multispectral
% filter-set proposals via band reduction + nested selection under LODO+k=3.
%
% KEY UPGRADES vs v1:
%   (1) Candidate pool derived directly from VIP profiles (local maxima with
%       VIP >= threshold), not restricted to "Top-5 windows".
%   (2) For FUSION, candidates are BLOCK-AWARE (FX10 vs FX17). Band simulation
%       is done within the corresponding block to avoid mixing in overlap.
%   (3) Candidate band responses are precomputed once per endpoint to speed
%       up greedy forward selection (SFS).
%
% INPUTS (must be in the current working directory or MATLAB path):
%   - Matriz_CHEM_HSI_MASTER_96.xlsx      (sheet 'Matriz')
%   - O5_VIP_Profiles_SelectedTraits.xlsx (sheets named as endpoint vars)
%
% OUTPUTS:
%   - Objetivo_6/Tables/O6_v2_CandidateBands_fromVIPPeaks.xlsx
%   - Objetivo_6/Tables/O6_v2_Multispectral_Performance_vsK.xlsx
%   - Objetivo_6/Tables/O6_v2_SelectedBandSets_Stability.xlsx
%   - Objetivo_6/Figures/O6_v2_<Endpoint>_<Mode>_Kcurve.*  (.fig + .png)
%   - Objetivo_6/Figures/Figure7_O6_Retention_TwoBaselines.* (.fig + .png)
%   - Objetivo_6/Logs/O6_v2_runlog.txt
%
% REQUIREMENTS:
%   - Statistics and Machine Learning Toolbox (plsregress, cvpartition)
% -------------------------------------------------------------------------

rng(20260226,'twister');

% =================== Configuration ===================
cfg.inMatrixXlsx      = 'Matriz_CHEM_HSI_MASTER_96.xlsx';
cfg.inMatrixSheet     = 'Matriz';
cfg.inVIPProfilesXlsx = 'O5_VIP_Profiles_SelectedTraits.xlsx';

cfg.outDir            = fullfile(pwd,'Objetivo_6');
cfg.figDir            = fullfile(cfg.outDir,'Figures');
cfg.tabDir            = fullfile(cfg.outDir,'Tables');
cfg.logDir            = fullfile(cfg.outDir,'Logs');

cfg.outerScheme       = 'LODO';          % 'LODO' or 'KFold'
cfg.outerGroupVar     = 'SamplingDate';  % leave-one-date-out unit
cfg.outerK            = 5;               % only used if outerScheme='KFold'
cfg.nRepeats          = 3;               % operational repeats (k=3)
cfg.rngSeedBase       = 20260226;

cfg.innerK            = 5;               % inner K-Fold for SFS scoring + LV choice
cfg.maxLV_ms          = 10;
cfg.maxLV_hsi         = 20;
cfg.useOneSE          = true;

cfg.K_list            = 3:12;            % number of bands to evaluate
cfg.FWHM_nm           = 20;              % assumed optical filter FWHM (nm)

% VIP-peak candidate pool
cfg.vipThresh         = 1.00;            % minimum VIP to consider a peak
cfg.maxCandPerBlock   = 25;              % cap per block
cfg.mergeTol_nm       = 10;              % merge peaks closer than this (within block)

% K recommendation rule (optional; used only in log)
cfg.maxRMSEIncreasePct = 10;             % choose minimal K within +10% RMSE of baseline

% Endpoints (must match sheet names in VIP profile file)
cfg.endpoints = { ...
    'TA_gL', ...
    'TSS_Brix', ...
    'TotalPhenolicPotential_mgkg', ...
    'TotalAnthocyanins_mgkg', ...
    'pH', ...
    'MalicAcid_gL' ...
    };

% =================== Figure style (Q1-friendly, BIG FONTS) ===================
cfg.fontName = 'Times New Roman';
cfg.fontAx   = 13;    % tick labels
cfg.fontLbl  = 16;    % axis labels
cfg.fontLeg  = 13;    % legend
cfg.fontAnn  = 12;    % annotations
cfg.figDPI   = 450;

% Wide figure for Figure 7
cfg.fig7Pos  = [100 100 1400 520];

cfg.bestRMSE_LODOk3 = struct();
cfg.bestRMSE_LODOk3.TA_gL = 0.642;
cfg.bestRMSE_LODOk3.TSS_Brix = 1.120;
cfg.bestRMSE_LODOk3.TotalPhenolicPotential_mgkg = 624.404;
cfg.bestRMSE_LODOk3.TotalAnthocyanins_mgkg = 40.996;
cfg.bestRMSE_LODOk3.MalicAcid_gL = 0.491;
cfg.bestRMSE_LODOk3.pH = 0.135;


% =================== Setup ===================
ensureDir(cfg.outDir); ensureDir(cfg.figDir); ensureDir(cfg.tabDir); ensureDir(cfg.logDir);
logFile = fullfile(cfg.logDir,'O6_v2_runlog.txt');
diary(logFile); diary on;

fprintf('--- O6 v2 start ---\n');
fprintf('Matrix: %s (sheet: %s)\nVIP profiles: %s\n\n', cfg.inMatrixXlsx, cfg.inMatrixSheet, cfg.inVIPProfilesXlsx);

% Load matrix
T = readtable(cfg.inMatrixXlsx, 'Sheet', cfg.inMatrixSheet);
mustHaveCols(T, [cfg.endpoints, {cfg.outerGroupVar}]);

% Extract spectral blocks
[FX10, FX17] = extractSpectralBlocks(T);

All_Cands   = table();
All_Results = table();
All_Stab    = table();

% =================== Main loop ===================
for e = 1:numel(cfg.endpoints)
    endpoint = cfg.endpoints{e};
    fprintf('\n====================================================\n');
    fprintf('Endpoint: %s\n', endpoint);

    % Read VIP profile sheet for this endpoint
    try
        VIP = readtable(cfg.inVIPProfilesXlsx, 'Sheet', endpoint);
    catch
        fprintf('WARNING: VIP profile sheet not found for %s. Skipping.\n', endpoint);
        continue;
    end
    mustHaveCols(VIP, {'wavelength_nm','VIP','Block'});

    % Infer mode from VIP profile blocks
    mode = inferModeFromVIPProfile(VIP);
    fprintf('Inferred mode: %s\n', mode);

    % Build candidate pool from VIP peaks (block-aware)
    cand = buildCandidatesFromVIPPeaks(VIP, cfg.vipThresh, cfg.maxCandPerBlock, cfg.mergeTol_nm);
    fprintf('Candidate pool size: %d (VIP>=%.2f; max %d per block)\n', height(cand), cfg.vipThresh, cfg.maxCandPerBlock);

    if height(cand) < min(cfg.K_list)
        fprintf('WARNING: Not enough candidates for K >= %d. Skipping.\n', min(cfg.K_list));
        continue;
    end

    % Precompute multispectral band responses (one column per candidate)
    F = precomputeCandidateFeatures(FX10, FX17, cand, cfg.FWHM_nm);

    % Baseline HSI RAW (full spectrum under mode)
    [X_full, ~] = buildX(FX10, FX17, mode);
    y = T.(endpoint);
    g = T.(cfg.outerGroupVar);

    fprintf('Computing RAW full-spectrum baseline under %s...\n', cfg.outerScheme);
    base = outerEvaluatePLSR(X_full, y, g, cfg, cfg.maxLV_hsi, 'BASELINE_HSI_RAW');
    fprintf('Baseline RMSE (HSI RAW): %.4g; R2: %.4f\n', base.RMSE, base.R2);

    % Evaluate K-band multispectral across K_list
    for K = cfg.K_list
        if K > height(cand); continue; end
        fprintf('  K=%d ...\n', K);

        ms = outerEvaluateMS(F, y, g, cand, K, cfg);

        % Store performance row
        resRow = table({endpoint}, {mode}, K, base.RMSE, ms.RMSE, ms.R2, ms.RMSE/base.RMSE, ...
            'VariableNames', {'EndpointVar','Mode','K','RMSE_HSI_RAW','RMSE_MS','R2_MS','RMSE_Ratio_MS_vs_HSI'});
        All_Results = [All_Results; resRow]; %#ok<AGROW>

        % Store stability summary
        stabRow = summariseBandStability(endpoint, mode, K, cand, ms.SelectedIdx_allRepeats);
        All_Stab = [All_Stab; stabRow]; %#ok<AGROW>
    end

    % Append candidates for this endpoint
    cand.EndpointVar = repmat({endpoint}, height(cand), 1);
    All_Cands = [All_Cands; cand]; %#ok<AGROW>

    % Plot K-curve for this endpoint (BIG FONTS)
    plotEndpointKCurve_v2(endpoint, mode, All_Results, base.RMSE, cfg);

    % Recommend minimal K within +x% baseline RMSE (or best achievable if none)
    recommendKForEndpoint_v2(endpoint, mode, All_Results, base.RMSE, cfg);
end

% =================== Summary Figure 7 (BIG FONTS) ===================
plotFigure7_TwoBaselines(All_Results, cfg);

% =================== Export ===================
writetable(All_Cands,   fullfile(cfg.tabDir, 'O6_v2_CandidateBands_fromVIPPeaks.xlsx'));
writetable(All_Results, fullfile(cfg.tabDir, 'O6_v2_Multispectral_Performance_vsK.xlsx'));
writetable(All_Stab,    fullfile(cfg.tabDir, 'O6_v2_SelectedBandSets_Stability.xlsx'));

fprintf('\n--- O6 v2 done ---\n');
diary off;

end

% ======================================================================
% Helper functions
% ======================================================================

function ensureDir(p)
if ~exist(p,'dir'); mkdir(p); end
end

function mustHaveCols(T, req)
vars = string(T.Properties.VariableNames);
req = string(req);
missing = req(~ismember(req, vars));
if ~isempty(missing)
    error('Missing required columns: %s', strjoin(cellstr(missing), ', '));
end
end

function [FX10, FX17] = extractSpectralBlocks(T)
vars = T.Properties.VariableNames;
fx10_cols = vars(startsWith(vars,'FX10_nm_'));
fx17_cols = vars(startsWith(vars,'FX17_nm_'));
if isempty(fx10_cols) || isempty(fx17_cols)
    error('Spectral columns not found. Expected prefixes FX10_nm_ and FX17_nm_.');
end
FX10.wl = parseWavelengths(fx10_cols);
FX17.wl = parseWavelengths(fx17_cols);

[FX10.wl, ix10] = sort(FX10.wl); FX10.cols = fx10_cols(ix10);
[FX17.wl, ix17] = sort(FX17.wl); FX17.cols = fx17_cols(ix17);

FX10.X = table2array(T(:, FX10.cols));
FX17.X = table2array(T(:, FX17.cols));
end

function wl = parseWavelengths(colNames)
wl = nan(numel(colNames),1);
for i = 1:numel(colNames)
    c = colNames{i};
    tok = regexp(c, 'nm_(\d+)p(\d+)$', 'tokens', 'once');
    if isempty(tok)
        error('Cannot parse wavelength from column name: %s', c);
    end
    wl(i) = str2double(tok{1}) + str2double(tok{2})/100;
end
end

function mode = inferModeFromVIPProfile(VIP)
blocks = unique(string(VIP.Block));
has10 = any(strcmpi(blocks,'FX10'));
has17 = any(strcmpi(blocks,'FX17'));
if has10 && has17
    mode = 'FUSION';
elseif has10
    mode = 'FX10';
elseif has17
    mode = 'FX17';
else
    mode = 'FX10';
end
end

function [X, wl] = buildX(FX10, FX17, mode)
switch upper(mode)
    case 'FX10'
        X = FX10.X; wl = FX10.wl;
    case 'FX17'
        X = FX17.X; wl = FX17.wl;
    case 'FUSION'
        X = [FX10.X, FX17.X];
        wl = [FX10.wl; FX17.wl];
    otherwise
        error('Unknown mode: %s', mode);
end
wl = wl(:);
end

function cand = buildCandidatesFromVIPPeaks(VIP, vipThresh, maxPerBlock, mergeTol)
blocks = unique(string(VIP.Block));
cand = table();

for b = 1:numel(blocks)
    blk = blocks(b);
    V = VIP(strcmpi(string(VIP.Block), blk), :);
    V = sortrows(V, 'wavelength_nm');

    wl = V.wavelength_nm(:);
    vv = V.VIP(:);

    n = numel(vv);
    if n < 3; continue; end

    idx = (2:n-1)';
    isPeak = vv(idx) >= vv(idx-1) & vv(idx) >= vv(idx+1) & vv(idx) >= vipThresh;
    peakIdx = idx(isPeak);
    if isempty(peakIdx); continue; end

    peakWl  = wl(peakIdx);
    peakVip = vv(peakIdx);

    [peakVip, ord] = sort(peakVip, 'desc');
    peakWl = peakWl(ord);

    keepWl  = [];
    keepVip = [];
    for i = 1:numel(peakWl)
        if isempty(keepWl)
            keepWl = peakWl(i); keepVip = peakVip(i);
        else
            if all(abs(keepWl - peakWl(i)) > mergeTol)
                keepWl(end+1,1) = peakWl(i); %#ok<AGROW>
                keepVip(end+1,1)= peakVip(i); %#ok<AGROW>
            end
        end
        if numel(keepWl) >= maxPerBlock; break; end
    end

    rankBlock = (1:numel(keepWl))';
    tmp = table(rankBlock, repmat(cellstr(blk), numel(keepWl), 1), keepWl, keepVip, ...
        'VariableNames', {'RankInBlock','Block','Peak_nm','Peak_VIP'});
    cand = [cand; tmp]; %#ok<AGROW>
end

[~,ord] = sort(cand.Peak_VIP, 'desc');
cand = cand(ord,:);
cand.RankGlobal = (1:height(cand))';
cand = movevars(cand, 'RankGlobal', 'Before', 1);
end

function F = precomputeCandidateFeatures(FX10, FX17, cand, FWHM)
half = FWHM/2;
n = size(FX10.X,1);
m = height(cand);
F = nan(n,m);

for j = 1:m
    blk = string(cand.Block(j));
    c   = cand.Peak_nm(j);

    if strcmpi(blk,'FX10')
        wl = FX10.wl; X = FX10.X;
    else
        wl = FX17.wl; X = FX17.X;
    end

    idx = find(wl >= (c-half) & wl <= (c+half));
    if isempty(idx)
        [~,k] = min(abs(wl - c));
        idx = k;
    end
    F(:,j) = mean(X(:,idx), 2);
end
end

function out = outerEvaluatePLSR(X, y, groupVar, cfg, maxLV, tag)
ok = ~isnan(y) & all(~isnan(X),2);
X = X(ok,:); y = y(ok); g = groupVar(ok);

n = size(X,1);
yHat_all = nan(n, cfg.nRepeats);

for r = 1:cfg.nRepeats
    rng(cfg.rngSeedBase + 1000*r, 'twister');
    yHat_all(:,r) = outerPredictPLSR_once(X, y, g, cfg, maxLV);
end

yHat = mean(yHat_all, 2, 'omitnan');
rmse = sqrt(mean((y - yHat).^2));
r2   = 1 - sum((y - yHat).^2) / sum((y - mean(y)).^2);

out = struct('Tag',tag,'RMSE',rmse,'R2',r2,'y',y,'yHat',yHat);
end

function yHat = outerPredictPLSR_once(X, y, g, cfg, maxLV)
n = size(X,1);
yHat = nan(n,1);
splits = makeOuterSplits(g, cfg);

for s = 1:numel(splits)
    te = splits{s}.testIdx;
    tr = splits{s}.trainIdx;

    Xtr = X(tr,:); ytr = y(tr);
    Xte = X(te,:);

    lv = selectLV_byInnerCV(Xtr, ytr, cfg.innerK, maxLV, cfg.useOneSE);

    [XtrZ, mu, sig] = zscoreFold(Xtr);
    XteZ = (Xte - mu) ./ sig;

    [~,~,~,~,beta] = plsregress(XtrZ, ytr, lv);
    yHat(te) = [ones(sum(te),1), XteZ] * beta;
end
end

function ms = outerEvaluateMS(F, y, groupVar, cand, K, cfg)
ok = ~isnan(y) & all(~isnan(F),2);
F = F(ok,:); y = y(ok); g = groupVar(ok);

n = size(F,1);
yHat_all = nan(n, cfg.nRepeats);
SelectedIdx_allRepeats = cell(cfg.nRepeats,1);

for r = 1:cfg.nRepeats
    rng(cfg.rngSeedBase + 2000*r, 'twister');
    [yHat_all(:,r), SelectedIdx_allRepeats{r}] = outerPredictMS_once(F, y, g, K, cfg);
end

yHat = mean(yHat_all, 2, 'omitnan');
rmse = sqrt(mean((y - yHat).^2));
r2   = 1 - sum((y - yHat).^2) / sum((y - mean(y)).^2);

ms = struct('RMSE',rmse,'R2',r2,'yHat',yHat,'SelectedIdx_allRepeats',{SelectedIdx_allRepeats});
end

function [yHat, selIdxCell] = outerPredictMS_once(F, y, g, K, cfg)
n = size(F,1);
yHat = nan(n,1);
splits = makeOuterSplits(g, cfg);
selIdxCell = cell(numel(splits),1);

for s = 1:numel(splits)
    te = splits{s}.testIdx;
    tr = splits{s}.trainIdx;

    Ftr = F(tr,:); ytr = y(tr);
    Fte = F(te,:);

    selIdx = greedySFS_features(Ftr, ytr, cfg, K);
    selIdxCell{s} = selIdx(:)';

    lv = selectLV_byInnerCV(Ftr(:,selIdx), ytr, cfg.innerK, cfg.maxLV_ms, cfg.useOneSE);

    [XtrZ, mu, sig] = zscoreFold(Ftr(:,selIdx));
    XteZ = (Fte(:,selIdx) - mu) ./ sig;

    [~,~,~,~,beta] = plsregress(XtrZ, ytr, lv);
    yHat(te) = [ones(sum(te),1), XteZ] * beta;
end
end

function selIdx = greedySFS_features(Ftr, ytr, cfg, K)
p = size(Ftr,2);
remaining = 1:p;
selected = [];

for step = 1:K
    bestIdx = NaN;
    bestRMSE = inf;

    for j = remaining
        trial = [selected, j];
        rmseCV = innerCV_RMSE_forSubset(Ftr(:,trial), ytr, cfg.innerK, cfg.maxLV_ms, cfg.useOneSE);
        if rmseCV < bestRMSE
            bestRMSE = rmseCV;
            bestIdx = j;
        end
    end

    selected(end+1) = bestIdx; %#ok<AGROW>
    remaining(remaining==bestIdx) = [];
end

selIdx = sort(selected(:))';
end

function rmseCV = innerCV_RMSE_forSubset(X, y, K, maxLV, useOneSE)
maxLV = min([maxLV, size(X,2), numel(y)-1]);
if maxLV < 1
    rmseCV = inf; return;
end

cv = cvpartition(numel(y), 'KFold', K);
rmseLV = nan(maxLV,1);
seLV   = nan(maxLV,1);

for lv = 1:maxLV
    foldRMSE = nan(cv.NumTestSets,1);
    for k = 1:cv.NumTestSets
        tr = training(cv,k); te = test(cv,k);
        Xtr = X(tr,:); ytr = y(tr);
        Xte = X(te,:); yte = y(te);

        [XtrZ, mu, sig] = zscoreFold(Xtr);
        XteZ = (Xte - mu) ./ sig;

        [~,~,~,~,beta] = plsregress(XtrZ, ytr, lv);
        yhat = [ones(sum(te),1), XteZ] * beta;
        foldRMSE(k) = sqrt(mean((yte - yhat).^2));
    end
    rmseLV(lv) = mean(foldRMSE,'omitnan');
    seLV(lv)   = std(foldRMSE,'omitnan')/sqrt(numel(foldRMSE));
end

[rmseMin, iMin] = min(rmseLV);

if ~useOneSE
    lvOpt = iMin;
else
    thresh = rmseMin + seLV(iMin);
    lvOpt = find(rmseLV <= thresh, 1, 'first');
    if isempty(lvOpt); lvOpt = iMin; end
end

rmseCV = rmseLV(lvOpt);
end

function splits = makeOuterSplits(g, cfg)
switch upper(cfg.outerScheme)
    case 'LODO'
        u = unique(g);
        splits = cell(numel(u),1);
        for i = 1:numel(u)
            te = (g == u(i));
            tr = ~te;
            splits{i} = struct('trainIdx',tr,'testIdx',te,'label',u(i));
        end
    case 'KFOLD'
        cv = cvpartition(numel(g), 'KFold', cfg.outerK);
        splits = cell(cv.NumTestSets,1);
        for i = 1:cv.NumTestSets
            te = test(cv,i);
            tr = training(cv,i);
            splits{i} = struct('trainIdx',tr,'testIdx',te,'label',i);
        end
    otherwise
        error('Unknown outerScheme: %s', cfg.outerScheme);
end
end

function lvOpt = selectLV_byInnerCV(X, y, K, maxLV, useOneSE)
maxLV = min([maxLV, size(X,2), numel(y)-1]);
if maxLV < 1
    lvOpt = 1; return;
end

cv = cvpartition(numel(y), 'KFold', K);
rmseLV = nan(maxLV,1);
seLV   = nan(maxLV,1);

for lv = 1:maxLV
    foldRMSE = nan(cv.NumTestSets,1);
    for k = 1:cv.NumTestSets
        tr = training(cv,k); te = test(cv,k);
        Xtr = X(tr,:); ytr = y(tr);
        Xte = X(te,:); yte = y(te);

        [XtrZ, mu, sig] = zscoreFold(Xtr);
        XteZ = (Xte - mu) ./ sig;

        [~,~,~,~,beta] = plsregress(XtrZ, ytr, lv);
        yhat = [ones(sum(te),1), XteZ] * beta;
        foldRMSE(k) = sqrt(mean((yte - yhat).^2));
    end
    rmseLV(lv) = mean(foldRMSE,'omitnan');
    seLV(lv)   = std(foldRMSE,'omitnan')/sqrt(numel(foldRMSE));
end

[rmseMin, iMin] = min(rmseLV);

if ~useOneSE
    lvOpt = iMin;
else
    thresh = rmseMin + seLV(iMin);
    lvOpt = find(rmseLV <= thresh, 1, 'first');
    if isempty(lvOpt); lvOpt = iMin; end
end
end

function [Xz, mu, sig] = zscoreFold(X)
mu = mean(X,1,'omitnan');
sig = std(X,0,1,'omitnan');
sig(sig==0) = 1;
Xz = (X - mu) ./ sig;
end

function stabRow = summariseBandStability(endpoint, mode, K, cand, SelectedIdx_allRepeats)
m = height(cand);
counts = zeros(m,1);

for r = 1:numel(SelectedIdx_allRepeats)
    selPerSplit = SelectedIdx_allRepeats{r};
    for s = 1:numel(selPerSplit)
        idx = selPerSplit{s};
        counts(idx) = counts(idx) + 1;
    end
end

[~,ord] = sort(counts,'descend');
ord = ord(1:min(K,numel(ord)));

bandStr = cell(numel(ord),1);
cntStr  = cell(numel(ord),1);
for i = 1:numel(ord)
    j = ord(i);
    bandStr{i} = sprintf('%s:%.2f', string(cand.Block(j)), cand.Peak_nm(j));
    cntStr{i}  = sprintf('%d', counts(j));
end

stabRow = table({endpoint}, {mode}, K, {strjoin(bandStr, ', ')}, {strjoin(cntStr, ', ')}, ...
    'VariableNames', {'EndpointVar','Mode','K','TopBands_ID','TopBands_counts'});
end

function plotEndpointKCurve_v2(endpoint, mode, All_Results, baseRMSE, cfg)
rows = All_Results(strcmp(All_Results.EndpointVar, endpoint) & strcmp(All_Results.Mode, mode), :);
if isempty(rows); return; end
rows = sortrows(rows, 'K');

fig = figure('Color','w');
ax = axes(fig); hold(ax,'on');

plot(ax, rows.K, rows.RMSE_MS, '-o', 'LineWidth', 1.5);
yline(ax, baseRMSE, '--', 'LineWidth', 1.2);

xlabel(ax, 'Number of bands (K)', 'FontName', cfg.fontName, 'FontSize', cfg.fontLbl);
ylabel(ax, 'RMSE', 'FontName', cfg.fontName, 'FontSize', cfg.fontLbl);
title(ax, sprintf('O6 v2: %s (%s)', endpoint, mode), 'FontName', cfg.fontName, 'FontSize', cfg.fontLbl);

set(ax, 'FontName', cfg.fontName, 'FontSize', cfg.fontAx, 'LineWidth', 1.1);
grid(ax,'on'); box(ax,'on');

ax.XLim = [min(rows.K) max(rows.K)];
setAxisTicksWithEnds(ax);

set(fig, 'Position', [100 100 1000 470]);

outBase = fullfile(cfg.figDir, sprintf('O6_v2_%s_%s_Kcurve', endpoint, mode));
savefig(fig, [outBase '.fig']);
exportgraphics(fig, [outBase '.png'], 'Resolution', cfg.figDPI);
close(fig);
end

function setAxisTicksWithEnds(ax)
xt = ax.XTick; xl = ax.XLim;
if isempty(xt); xt = linspace(xl(1), xl(2), 5); end
xt = unique([xl(1), xt, xl(2)]);
ax.XTick = xt;

yt = ax.YTick; yl = ax.YLim;
if isempty(yt); yt = linspace(yl(1), yl(2), 5); end
yt = unique([yl(1), yt, yl(2)]);
ax.YTick = yt;
end

function recommendKForEndpoint_v2(endpoint, mode, All_Results, baseRMSE, cfg)
rows = All_Results(strcmp(All_Results.EndpointVar, endpoint) & strcmp(All_Results.Mode, mode), :);
if isempty(rows); return; end
rows = sortrows(rows, 'K');

thresh = baseRMSE * (1 + cfg.maxRMSEIncreasePct/100);
ok = rows.RMSE_MS <= thresh;

if any(ok)
    Krec = rows.K(find(ok,1,'first'));
    fprintf('Recommendation: K=%d (<= +%.1f%% baseline RMSE)\n', Krec, cfg.maxRMSEIncreasePct);
else
    [~,iBest] = min(rows.RMSE_MS);
    fprintf('Recommendation: none within +%.1f%% baseline. Best RMSE at K=%d.\n', cfg.maxRMSEIncreasePct, rows.K(iBest));
end
end

% =================== FIGURE 7 (BIG FONTS) ===================

function plotFigure7_TwoBaselines(All_Results, cfg)
if isempty(All_Results) || height(All_Results)==0
    fprintf('Figure 7: All_Results is empty.\n');
    return;
end

endpoints = unique(All_Results.EndpointVar, 'stable');
S = table();

for i = 1:numel(endpoints)
    ep = endpoints{i};
    rows = All_Results(strcmp(All_Results.EndpointVar, ep), :);
    if isempty(rows); continue; end

    [~,ib] = min(rows.RMSE_MS);
    r = rows(ib,:);

    ratioRAW = r.RMSE_MS / r.RMSE_HSI_RAW;

    ratioBEST = NaN;
    if isfield(cfg,'bestRMSE_LODOk3') && isstruct(cfg.bestRMSE_LODOk3) && isfield(cfg.bestRMSE_LODOk3, ep)
        baseBest = cfg.bestRMSE_LODOk3.(ep);
        if ~isempty(baseBest) && ~isnan(baseBest) && baseBest > 0
            ratioBEST = r.RMSE_MS / baseBest;
        end
    end

    S = [S; table({ep}, r.Mode, r.K, ratioRAW, ratioBEST, ...
        'VariableNames', {'EndpointVar','Mode','K','RatioRAW','RatioBEST'})]; %#ok<AGROW>
end

if isempty(S)
    fprintf('Figure 7: no endpoints aggregated.\n');
    return;
end

labels = cell(height(S),1);
for i = 1:height(S)
    labels{i} = prettyEndpointLabel(S.EndpointVar{i});
end

fig = figure('Color','w');
set(fig,'Units','pixels','Position',cfg.fig7Pos);
ax = axes(fig); hold(ax,'on');

x = 1:height(S);
w = 0.36;

bar(ax, x - w/2, S.RatioRAW, w, 'DisplayName','vs full-spectrum RAW');

hasBEST = any(~isnan(S.RatioBEST));
if hasBEST
    bar(ax, x + w/2, S.RatioBEST, w, 'DisplayName','vs best full-spectrum (LODO, k=3)');
end

yline(ax, 1.0, '--', 'LineWidth', 1.2);

set(ax,'FontName',cfg.fontName,'FontSize',cfg.fontAx,'LineWidth',1.1);
ax.Box = 'on';

ax.XTick = x;
ax.XTickLabel = labels;

ylabel(ax,'RMSE ratio (multispectral / baseline)','FontName',cfg.fontName,'FontSize',cfg.fontLbl);
xlabel(ax,'Endpoint','FontName',cfg.fontName,'FontSize',cfg.fontLbl);

lgd = legend(ax,'Location','northwest');
lgd.FontSize = 11;          % baja desde 13 si no cabe
lgd.NumColumns = 1;         % o 2 si te interesa compacto horizontal
lgd.Box = 'on';

for i = 1:height(S)
    ytop = S.RatioRAW(i);
    if hasBEST && ~isnan(S.RatioBEST(i))
        ytop = max(ytop, S.RatioBEST(i));
    end
    text(ax, x(i), ytop + 0.06, sprintf('K=%d', S.K(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontName',cfg.fontName,'FontSize',cfg.fontAnn);
end

ymax = max([S.RatioRAW; S.RatioBEST], [], 'omitnan');
ax.YLim = [0, max(1.25, ymax*1.12)];

outBase = fullfile(cfg.figDir, 'Figure7_O6_Retention_TwoBaselines');
savefig(fig, [outBase '.fig']);
exportgraphics(fig, [outBase '.png'], 'Resolution', cfg.figDPI);
close(fig);

fprintf('Figure 7 saved to: %s(.fig/.png)\n', outBase);
end

function s = prettyEndpointLabel(ep)
switch ep
   case 'TA_gL', s = 'TA';
case 'TSS_Brix', s = 'TSS';
case 'TotalPhenolicPotential_mgkg', s = 'PPT';
case 'TotalAnthocyanins_mgkg', s = 'AnT';
case 'MalicAcid_gL', s = 'Malic';
case 'pH', s = 'pH';
    otherwise
        s = ep;
end
end