function O4_LODO_Strict_and_k3_run_v2()
% O4_LODO_Strict_and_k3_run_v2
% ------------------------------------------------------------
% Objective 4 / Results 3.4:
% LODO (Leave-One-Date-Out; group = SamplingEvent) evaluated under:
%   (A) Strict LODO (no information from the held-out date)
%   (B) LODO + k anchors/date (minimal daily recalibration)
%
% Anchor strategy (deployment-realistic):
%   - For each held-out date, select k "anchor" samples using spectral diversity
%     (Kennard–Stone on preprocessed Xtest).
%   - Assume these k anchors are chemically assayed on that date.
%   - Compute intercept update: bias = mean(y_true_anchor - y_pred_anchor)
%   - Apply y_pred_adj = y_pred + bias to the remaining samples of that date.
%   - Evaluate performance ONLY on non-anchor test samples (to avoid leakage).
%
% INPUT:
%   - Matriz_CHEM_HSI_MASTER_96.xlsx (Sheets: 'Matriz', 'DataDictionary')
%
% OUTPUT (./Objetivo_4/):
%   - LODO_Performance_SelectedModels_Strict_vs_kAnchors.xlsx
%       * Sheet 'Table5' : compact summary for manuscript Table 5
%       * Sheet 'PerFold': per-held-out-date diagnostics
%   - LODO_Predictions_SelectedModels_Strict_vs_kAnchors.xlsx
%       * Sheet 'Predictions': long-format predictions with AnchorFlag
%   - Figure6_LODO_Strict_vs_kAnchors.{fig,png} : slope chart (R² strict vs k-anchors)
%
% Notes:
%   - Selected configurations (sensor × preprocessing × LV, plus R²cv_random)
%     are hard-coded from Table 4 (as used in your pipeline).
%   - No outlier filtering is applied here (kept minimal and deterministic).
%
% ------------------------------------------------------------

rng(20260225,'twister');

% =================== Configuration ===================
cfg = struct();
cfg.inputFile  = fullfile(pwd,'Matriz_CHEM_HSI_MASTER_96.xlsx');
cfg.sheetData  = 'Matriz';
cfg.sheetDict  = 'DataDictionary';

cfg.outDir     = fullfile(pwd,'Objetivo_4');
cfg.figureDPI  = 400;

% Anchor calibration
cfg.kAnchors   = 3;     % <-- the "k=3" you requested
cfg.anchorMode = 'KS';  % 'KS' (Kennard–Stone) or 'FirstK' (by LabCode order)

% SG settings (match your pipeline defaults)
cfg.sgPolyOrder = 2;
cfg.sgFrameLen  = 15;

% Fusion rule
cfg.fusionCut_nm = 950; % FX10 <= 950, FX17 > 950
cfg.fusionL2norm = true;

% =================== Read data ===================
assert(isfile(cfg.inputFile), 'Input file not found: %s', cfg.inputFile);
T  = readtable(cfg.inputFile, 'Sheet', cfg.sheetData, 'PreserveVariableNames', true);
DD = readtable(cfg.inputFile, 'Sheet', cfg.sheetDict, 'PreserveVariableNames', true);

% Group variable (LODO)
if ismember('SamplingEvent', T.Properties.VariableNames)
    groupVar = 'SamplingEvent';
elseif ismember('SamplingDate', T.Properties.VariableNames)
    groupVar = 'SamplingDate';
else
    error('Neither SamplingEvent nor SamplingDate found in the matrix.');
end

% Identify spectral blocks
fx10Cols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'FX10_nm_'));
fx17Cols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'FX17_nm_'));
assert(~isempty(fx10Cols) && ~isempty(fx17Cols), 'FX10/FX17 spectral columns not found.');

fx10_nm = parseNmFromVarNames(fx10Cols, 'FX10_nm_');
fx17_nm = parseNmFromVarNames(fx17Cols, 'FX17_nm_');

% =================== Selected models from Table 4 ===================
sel = table();

% Endpoint variables must match the Excel matrix column names:
sel.EndpointVar = { ...
    'BerryWeight50_g'
    'TSS_Brix'
    'pH'
    'TA_gL'
    'MalicAcid_gL'
    'TartaricAcid_gL'
    'Potassium_gL'
    'AlphaAminoN_mgL'
    'Ammonium_mgL'
    'TotalPhenolicPotential_mgkg'
    'ExtractablePhenolicPotential_mgkg'
    'TotalAnthocyanins_mgkg'
    'ExtractableAnthocyanins_mgkg'
    };
sel.EndpointVar = sel.EndpointVar(:);

sel.Sensor = { ...
    'FUSION'
    'FUSION'
    'FUSION'
    'FX10'
    'FUSION'
    'FUSION'
    'FX17'
    'FX17'
    'FUSION'
    'FX10'
    'FX17'
    'FUSION'
    'FUSION'
    };
sel.Sensor = sel.Sensor(:);

sel.Preproc = { ...
    'SNV+SG1'
    'RAW'
    'SG1'
    'SG2'
    'SNV'
    'RAW'
    'SG2'
    'SNV'
    'SG2'
    'SG2'
    'SG1'
    'SG2'
    'SG1'
    };
sel.Preproc = sel.Preproc(:);

% IMPORTANT: must be a COLUMN vector (13x1)
sel.LV = [ ...
     3
     2
     3
     5
     5
     2
     2
     1
     8
     7
     1
     2
     2
    ];
sel.LV = sel.LV(:);

% Random CV baseline from Table 4 (must be COLUMN vector)
sel.R2cv_random = [ ...
    0.57
    0.76
    0.67
    0.85
    0.70
    0.40
    0.16
    0.07
    0.21
    0.79
    0.46
    0.64
    0.08
    ];
sel.R2cv_random = sel.R2cv_random(:);

assert(height(sel)==numel(sel.LV) && height(sel)==numel(sel.R2cv_random), 'Selection table length mismatch.');

% Add human-readable endpoint label (with units) from DataDictionary
sel = addUnitsAndLabels(sel, DD);

% =================== Main loop ===================
groups = unique(T.(groupVar));
nG = numel(groups);

perfRows = [];
perFoldRows = [];
predRows = [];

for i = 1:height(sel)

    yVar   = sel.EndpointVar{i};

    sensor = sel.Sensor{i};

    prep   = sel.Preproc{i};

    LV     = sel.LV(i);

    assert(ismember(yVar, T.Properties.VariableNames), 'Endpoint not found: %s', yVar);

    yAll = T.(yVar);
    okY  = isfinite(yAll);

    % Storage for pooled evaluation across folds:
    yTrue_strict = [];
    yPred_strict = [];
    yTrue_k      = [];
    yPred_k      = [];

    for g = 1:nG
        gVal = groups(g);

        isTest  = okY & ismember(T.(groupVar), gVal);
        isTrain = okY & ~ismember(T.(groupVar), gVal);

        if ~any(isTest) || ~any(isTrain)
            continue;
        end

        % Build X blocks
        [Xtrain, Xtest] = buildXsensor(T, isTrain, isTest, sensor, ...
            fx10Cols, fx17Cols, fx10_nm, fx17_nm, cfg);

        yTrain = yAll(isTrain);
        yTest  = yAll(isTest);

        % Preprocess (deterministic)
        [XtrainP, XtestP] = applyPreproc(Xtrain, Xtest, sensor, prep, ...
            fx10_nm, fx17_nm, cfg);

        % Fit final PLS-R (fixed LV)
        LVfit = min([LV, size(XtrainP,1)-1, rank(XtrainP)]);
        if LVfit < 1
            warning('Endpoint %s: LVfit < 1 for fold %d. Skipping fold.', yVar, g);
            continue;
        end

        [~,~,~,~,beta] = plsregress(XtrainP, yTrain, LVfit);
        yHat = [ones(size(XtestP,1),1) XtestP] * beta;

        % ---------- Strict LODO metrics (on ALL test samples) ----------
        [r2_s, rmse_s, bias_s] = regressionMetrics(yTest, yHat);

        % Pool for strict
        yTrue_strict = [yTrue_strict; yTest(:)]; %#ok<AGROW>
        yPred_strict = [yPred_strict; yHat(:)];  %#ok<AGROW>

        % ---------- k-anchors/date intercept update ----------
        k = cfg.kAnchors;
        nTest = numel(yTest);

        % Select anchors based on XtestP only (deployment-realistic)
        anchorIdx = selectAnchors(XtestP, T, isTest, k, cfg);

        % Compute bias on anchors and adjust others
        yHat_adj = yHat;
        if ~isempty(anchorIdx)
            bias_k = mean(yTest(anchorIdx) - yHat(anchorIdx), 'omitnan');
            yHat_adj = yHat + bias_k;
        else
            bias_k = NaN;
        end

        % Evaluate ONLY on non-anchor samples
        evalMask = true(nTest,1);
        evalMask(anchorIdx) = false;

        [r2_k, rmse_k, bias_k_eval] = regressionMetrics(yTest(evalMask), yHat_adj(evalMask));

        % Pool for k-anchors evaluation (non-anchors only)
        yTrue_k = [yTrue_k; yTest(evalMask)]; %#ok<AGROW>
        yPred_k = [yPred_k; yHat_adj(evalMask)]; %#ok<AGROW>

        % Per-fold diagnostics
        perFoldRows = [perFoldRows; { ...
            sel.Endpoint{i}, yVar, sensor, prep, LVfit, groupToString(gVal), ...
            nTest, ...
            r2_s, rmse_s, bias_s, ...
            r2_k, rmse_k, bias_k_eval, ...
            numel(anchorIdx), bias_k ...
            }]; %#ok<AGROW>

        % Prediction rows (long-format) — include both strict and adjusted
        % AnchorFlag indicates samples used for daily recalibration.
        idxTest = find(isTest);
        for j = 1:numel(idxTest)
            r = idxTest(j);
            isAnchor = ismember(j, anchorIdx);
            predRows = [predRows; { ...
                sel.Endpoint{i}, yVar, sensor, prep, LV, ...
                labCodeToString(T, r), ...
                groupToString(T.(groupVar)(r)), ...
                dateToString(T, r), ...
                string(T.IrrigationRegime(r)), ...
                string(T.VineyardFloorManagement(r)), ...
                string(T.Block(r)), ...
                double(yAll(r)), ...
                double(yHat(j)), ...
                double(yHat_adj(j)), ...
                double(yHat(j)-yAll(r)), ...
                double(yHat_adj(j)-yAll(r)), ...
                logical(isAnchor) ...
                }]; %#ok<AGROW>
        end

    end % folds

    % Pooled metrics
    [R2_LODO_strict, RMSE_LODO_strict, ~] = regressionMetrics(yTrue_strict, yPred_strict);
    [R2_LODO_k,      RMSE_LODO_k,      ~] = regressionMetrics(yTrue_k, yPred_k);

    dR2_strict = R2_LODO_strict - sel.R2cv_random(i);
    dR2_k      = R2_LODO_k      - sel.R2cv_random(i);

    perfRows = [perfRows; { ...
        sel.Endpoint{i}, ...
        yVar, ...
        sprintf('%s-%s (LV=%d)', sensor, prep, LV), ...
        sel.R2cv_random(i), ...
        R2_LODO_strict, dR2_strict, RMSE_LODO_strict, ...
        R2_LODO_k,      dR2_k,      RMSE_LODO_k, ...
        sum(okY) ...
        }]; %#ok<AGROW>

end % endpoints

% =================== Build output tables ===================
perfTable = cell2table(perfRows, 'VariableNames', { ...
    'Endpoint', 'EndpointVar', 'BestConfig', 'R2cv_random', ...
    'R2_LODO_strict', 'DeltaR2_strict', 'RMSE_LODO_strict', ...
    sprintf('R2_LODO_k%d', cfg.kAnchors), sprintf('DeltaR2_k%d', cfg.kAnchors), sprintf('RMSE_LODO_k%d', cfg.kAnchors), ...
    'n_total'});

perFoldTable = cell2table(perFoldRows, 'VariableNames', { ...
    'Endpoint','EndpointVar','Sensor','Preproc','LV_fit','HeldOutGroup','nTest', ...
    'R2_strict','RMSE_strict','Bias_strict', ...
    sprintf('R2_k%d', cfg.kAnchors), sprintf('RMSE_k%d', cfg.kAnchors), sprintf('Bias_k%d_eval', cfg.kAnchors), ...
    'kAnchors','Bias_interceptUpdate'});

predTable = cell2table(predRows, 'VariableNames', { ...
    'Endpoint','EndpointVar','Sensor','Preproc','LV_Table4', ...
    'LabCode','HeldOutGroup','SamplingDate','IrrigationRegime','VineyardFloorManagement','Block', ...
    'Y_true','Y_pred_strict','Y_pred_k','Residual_strict','Residual_k','AnchorFlag'});

% =================== Write outputs ===================
if ~exist(cfg.outDir,'dir'); mkdir(cfg.outDir); end

outPerf = fullfile(cfg.outDir, sprintf('LODO_Performance_SelectedModels_Strict_vs_k%dAnchors.xlsx', cfg.kAnchors));
outPred = fullfile(cfg.outDir, sprintf('LODO_Predictions_SelectedModels_Strict_vs_k%dAnchors.xlsx', cfg.kAnchors));

writetable(perfTable, outPerf, 'Sheet', 'Table5');
writetable(perFoldTable, outPerf, 'Sheet', 'PerFold');
writetable(predTable, outPred, 'Sheet', 'Predictions');

% =================== Figure 6: slope chart (R² strict vs R² k) ===================
fig = figure('Color','w','Position',[100 100 1200 520]);
tPlot = perfTable;
% Sort by k-anchors R² (ascending for readability)
kName = sprintf('R2_LODO_k%d', cfg.kAnchors);
tPlot = sortrows(tPlot, kName, 'ascend');

y = (1:height(tPlot))';

x1 = tPlot.R2_LODO_strict;
x2 = tPlot.(kName);

plot(x1, y, 'o-', 'LineWidth', 1.2); hold on;
plot(x2, y, 'o-', 'LineWidth', 1.2);

% Connectors
for ii=1:height(tPlot)
    plot([x1(ii) x2(ii)], [y(ii) y(ii)], '-', 'LineWidth', 0.8);
end

grid on;
xline(0,'-');
set(gca,'YTick',y,'YTickLabel',tPlot.Endpoint,'FontName','Times New Roman','FontSize',9);
xlabel('Pooled R^2','FontName','Times New Roman','FontSize',12);
ylabel('Endpoint','FontName','Times New Roman','FontSize',12);
title(sprintf('Objective 4 – LODO robustness (strict) vs LODO + k=%d anchors/date', cfg.kAnchors), ...
    'FontName','Times New Roman','FontSize',12);

legend({'Strict LODO (k=0)','LODO + anchors (intercept update)'}, 'Location','best', 'FontSize',9);

figFig = fullfile(cfg.outDir, sprintf('Figure6_LODO_Strict_vs_k%dAnchors.fig', cfg.kAnchors));
figPng = fullfile(cfg.outDir, sprintf('Figure6_LODO_Strict_vs_k%dAnchors.png', cfg.kAnchors));
savefig(fig, figFig);
exportgraphics(fig, figPng, 'Resolution', cfg.figureDPI);

% =================== Run log ===================
logFile = fullfile(cfg.outDir, sprintf('O4_LODO_k%d_runlog.txt', cfg.kAnchors));
fid = fopen(logFile,'w');
fprintf(fid, 'O4 LODO robustness: strict vs k-anchors/date\n');
fprintf(fid, 'Timestamp: %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'Input: %s\n', cfg.inputFile);
fprintf(fid, 'Group variable: %s\n', groupVar);
fprintf(fid, 'kAnchors: %d (mode=%s)\n', cfg.kAnchors, cfg.anchorMode);
fprintf(fid, 'SG: poly=%d, frame=%d\n', cfg.sgPolyOrder, cfg.sgFrameLen);
fprintf(fid, 'Fusion cut: FX10 <= %.1f nm; FX17 > %.1f nm; L2norm=%d\n', cfg.fusionCut_nm, cfg.fusionCut_nm, cfg.fusionL2norm);
fprintf(fid, 'Outputs:\n  %s\n  %s\n  %s\n', outPerf, outPred, figPng);
fclose(fid);

disp('--- Objective 4 complete ---');
disp(['Saved: ' outPerf]);
disp(['Saved: ' outPred]);
disp(['Saved: ' figPng]);

end

% =================== Helper functions ===================

function sel = addUnitsAndLabels(sel, DD)
unit = strings(height(sel),1);
desc = strings(height(sel),1);

if ~all(ismember({'Variable','Unit (recommended)','Description (EN)'}, DD.Properties.VariableNames))
    % Fall back (do not fail)
    sel.Endpoint = sel.EndpointVar;
    return;
end

for i=1:height(sel)
    v = sel.EndpointVar{i};
    hit = strcmp(DD.Variable, v);
    if any(hit)
        u = DD.("Unit (recommended)")(find(hit,1,'first'));
        d = DD.("Description (EN)")(find(hit,1,'first'));
        if ismissing(u); u = ""; end
        if ismissing(d); d = v; end
        unit(i) = string(u);
        desc(i) = string(d);
    else
        unit(i) = "";
        desc(i) = v;
    end
end

endpointLabel = desc;
for i=1:numel(endpointLabel)
    if strlength(unit(i))>0
        endpointLabel(i) = endpointLabel(i) + " (" + unit(i) + ")";
    end
end

sel.Endpoint = cellstr(endpointLabel);
end

function nm = parseNmFromVarNames(varNames, prefix)
nm = nan(numel(varNames),1);
for i=1:numel(varNames)
    s = varNames{i};
    s = erase(s, prefix);
    s = strrep(s,'p','.');
    nm(i) = str2double(s);
end
end

function [Xtrain, Xtest] = buildXsensor(T, isTrain, isTest, sensor, fx10Cols, fx17Cols, fx10_nm, fx17_nm, cfg)
switch upper(sensor)
    case 'FX10'
        Xtrain = T{isTrain, fx10Cols};
        Xtest  = T{isTest,  fx10Cols};

    case 'FX17'
        Xtrain = T{isTrain, fx17Cols};
        Xtest  = T{isTest,  fx17Cols};

    case 'FUSION'
        fx10Keep = fx10_nm <= cfg.fusionCut_nm;
        fx17Keep = fx17_nm >  cfg.fusionCut_nm;

        X10_tr = T{isTrain, fx10Cols(fx10Keep)};
        X10_te = T{isTest,  fx10Cols(fx10Keep)};
        X17_tr = T{isTrain, fx17Cols(fx17Keep)};
        X17_te = T{isTest,  fx17Cols(fx17Keep)};

        Xtrain = [X10_tr, X17_tr];
        Xtest  = [X10_te, X17_te];

    otherwise
        error('Unknown sensor: %s', sensor);
end
end

function Xn = l2normRows(X)
den = sqrt(sum(X.^2,2));
den(den==0) = eps;
Xn = X ./ den;
end

function [XtrainP, XtestP] = applyPreproc(Xtrain, Xtest, sensor, prep, fx10_nm, fx17_nm, cfg)
prep = upper(strrep(prep,' ','')); 
sensor = upper(sensor);

if strcmp(sensor,'FUSION')
    fx10Keep = fx10_nm <= cfg.fusionCut_nm;
    fx17Keep = fx17_nm >  cfg.fusionCut_nm;
    n10 = sum(fx10Keep);
    n17 = sum(fx17Keep);

    X10_tr = Xtrain(:,1:n10);
    X17_tr = Xtrain(:,n10+1:n10+n17);
    X10_te = Xtest(:,1:n10);
    X17_te = Xtest(:,n10+1:n10+n17);

    d10 = median(diff(fx10_nm(fx10Keep)), 'omitnan');
    d17 = median(diff(fx17_nm(fx17Keep)), 'omitnan');

    switch prep
        case 'RAW'
        case 'SNV'
            X10_tr = snvRows(X10_tr);  X10_te = snvRows(X10_te);
            X17_tr = snvRows(X17_tr);  X17_te = snvRows(X17_te);
        case 'SG1'
            X10_tr = sgolayDerivative2D(X10_tr, cfg.sgPolyOrder, cfg.sgFrameLen, 1, d10);
            X10_te = sgolayDerivative2D(X10_te, cfg.sgPolyOrder, cfg.sgFrameLen, 1, d10);
            X17_tr = sgolayDerivative2D(X17_tr, cfg.sgPolyOrder, cfg.sgFrameLen, 1, d17);
            X17_te = sgolayDerivative2D(X17_te, cfg.sgPolyOrder, cfg.sgFrameLen, 1, d17);
        case 'SG2'
            X10_tr = sgolayDerivative2D(X10_tr, cfg.sgPolyOrder, cfg.sgFrameLen, 2, d10);
            X10_te = sgolayDerivative2D(X10_te, cfg.sgPolyOrder, cfg.sgFrameLen, 2, d10);
            X17_tr = sgolayDerivative2D(X17_tr, cfg.sgPolyOrder, cfg.sgFrameLen, 2, d17);
            X17_te = sgolayDerivative2D(X17_te, cfg.sgPolyOrder, cfg.sgFrameLen, 2, d17);
        case 'SNV+SG1'
            X10_tr = sgolayDerivative2D(snvRows(X10_tr), cfg.sgPolyOrder, cfg.sgFrameLen, 1, d10);
            X10_te = sgolayDerivative2D(snvRows(X10_te), cfg.sgPolyOrder, cfg.sgFrameLen, 1, d10);
            X17_tr = sgolayDerivative2D(snvRows(X17_tr), cfg.sgPolyOrder, cfg.sgFrameLen, 1, d17);
            X17_te = sgolayDerivative2D(snvRows(X17_te), cfg.sgPolyOrder, cfg.sgFrameLen, 1, d17);
        case 'SNV+SG2'
            X10_tr = sgolayDerivative2D(snvRows(X10_tr), cfg.sgPolyOrder, cfg.sgFrameLen, 2, d10);
            X10_te = sgolayDerivative2D(snvRows(X10_te), cfg.sgPolyOrder, cfg.sgFrameLen, 2, d10);
            X17_tr = sgolayDerivative2D(snvRows(X17_tr), cfg.sgPolyOrder, cfg.sgFrameLen, 2, d17);
            X17_te = sgolayDerivative2D(snvRows(X17_te), cfg.sgPolyOrder, cfg.sgFrameLen, 2, d17);
        otherwise
            error('Unknown preprocessing: %s', prep);
    end

    if cfg.fusionL2norm
        X10_tr = l2normRows(X10_tr);  X10_te = l2normRows(X10_te);
        X17_tr = l2normRows(X17_tr);  X17_te = l2normRows(X17_te);
    end

    XtrainP = [X10_tr, X17_tr];
    XtestP  = [X10_te, X17_te];
    return;
end

switch prep
    case 'RAW'
        XtrainP = Xtrain;
        XtestP  = Xtest;
    case 'SNV'
        XtrainP = snvRows(Xtrain);
        XtestP  = snvRows(Xtest);
    case 'SG1'
        if strcmp(sensor,'FX10'); delta = median(diff(fx10_nm), 'omitnan'); else; delta = median(diff(fx17_nm), 'omitnan'); end
        XtrainP = sgolayDerivative2D(Xtrain, cfg.sgPolyOrder, cfg.sgFrameLen, 1, delta);
        XtestP  = sgolayDerivative2D(Xtest,  cfg.sgPolyOrder, cfg.sgFrameLen, 1, delta);
    case 'SG2'
        if strcmp(sensor,'FX10'); delta = median(diff(fx10_nm), 'omitnan'); else; delta = median(diff(fx17_nm), 'omitnan'); end
        XtrainP = sgolayDerivative2D(Xtrain, cfg.sgPolyOrder, cfg.sgFrameLen, 2, delta);
        XtestP  = sgolayDerivative2D(Xtest,  cfg.sgPolyOrder, cfg.sgFrameLen, 2, delta);
    case 'SNV+SG1'
        if strcmp(sensor,'FX10'); delta = median(diff(fx10_nm), 'omitnan'); else; delta = median(diff(fx17_nm), 'omitnan'); end
        XtrainP = sgolayDerivative2D(snvRows(Xtrain), cfg.sgPolyOrder, cfg.sgFrameLen, 1, delta);
        XtestP  = sgolayDerivative2D(snvRows(Xtest),  cfg.sgPolyOrder, cfg.sgFrameLen, 1, delta);
    case 'SNV+SG2'
        if strcmp(sensor,'FX10'); delta = median(diff(fx10_nm), 'omitnan'); else; delta = median(diff(fx17_nm), 'omitnan'); end
        XtrainP = sgolayDerivative2D(snvRows(Xtrain), cfg.sgPolyOrder, cfg.sgFrameLen, 2, delta);
        XtestP  = sgolayDerivative2D(snvRows(Xtest),  cfg.sgPolyOrder, cfg.sgFrameLen, 2, delta);
    otherwise
        error('Unknown preprocessing: %s', prep);
end
end

function X = snvRows(X)
mu = mean(X,2,'omitnan');
sd = std(X,0,2,'omitnan');
sd(sd==0) = eps;
X = (X - mu) ./ sd;
end

function Xd = sgolayDerivative2D(X, polyOrder, frameLen, derivOrder, delta)
if mod(frameLen,2) ~= 1
    error('SG frameLen must be odd.');
end
halfwin = (frameLen-1)/2;
[~, g] = sgolay(polyOrder, frameLen);
filt = fliplr(g(:, derivOrder+1)'); % correlation-equivalent
left  = X(:, halfwin+1:-1:2);
right = X(:, end-1:-1:end-halfwin);
Xpad  = [left, X, right];
Xd = conv2(Xpad, filt, 'valid');
Xd = Xd * factorial(derivOrder) / (delta^derivOrder);
end

function anchorIdx = selectAnchors(XtestP, T, isTest, k, cfg)
n = size(XtestP,1);
k = min(k, n);
if k <= 0
    anchorIdx = [];
    return;
end

switch upper(cfg.anchorMode)
    case 'FIRSTK'
        anchorIdx = (1:k)';
    case 'KS'
        % Kennard–Stone on z-scored features (within test set)
        X = XtestP;
        mu = mean(X,1,'omitnan');
        sd = std(X,0,1,'omitnan');
        sd(sd==0) = 1;
        Z = (X - mu) ./ sd;
        anchorIdx = kennardStone(Z, k);
    otherwise
        anchorIdx = (1:k)';
end
anchorIdx = anchorIdx(:);
end

function idxSel = kennardStone(X, k)
% Simple Kennard–Stone selection (greedy farthest point)
% X: n x p
n = size(X,1);
if k >= n
    idxSel = (1:n)';
    return;
end

% Pairwise distances (n is small: ~24)
D = squareform(pdist(X,'euclidean'));

% Start with the farthest pair
[~, ind] = max(D(:));
[i1, i2] = ind2sub([n n], ind);
idxSel = [i1; i2];

while numel(idxSel) < k
    remaining = setdiff(1:n, idxSel);
    % For each remaining point, compute distance to nearest selected point
    minDist = zeros(numel(remaining),1);
    for r = 1:numel(remaining)
        rr = remaining(r);
        minDist(r) = min(D(rr, idxSel));
    end
    [~, j] = max(minDist);
    idxSel = [idxSel; remaining(j)]; %#ok<AGROW>
end
idxSel = idxSel(1:k);
idxSel = idxSel(:);
end

function [R2, RMSE, Bias] = regressionMetrics(yTrue, yPred)
yTrue = yTrue(:);
yPred = yPred(:);
idx = isfinite(yTrue) & isfinite(yPred);
yTrue = yTrue(idx);
yPred = yPred(idx);

if numel(yTrue) < 3
    R2 = NaN; RMSE = NaN; Bias = NaN; return;
end

res = yPred - yTrue;
SSE = sum(res.^2);
SST = sum((yTrue - mean(yTrue)).^2);

if SST <= 0
    R2 = NaN;
else
    R2 = 1 - SSE/SST;
end

RMSE = sqrt(mean(res.^2));
Bias = mean(res);
end

function s = groupToString(gVal)
if isdatetime(gVal)
    s = datestr(gVal,'yyyy-mm-dd');
else
    s = string(gVal);
end
s = char(s);
end

function s = dateToString(T, r)
if ismember('SamplingDate', T.Properties.VariableNames)
    d = T.SamplingDate(r);
    if isdatetime(d)
        s = datestr(d,'yyyy-mm-dd');
    else
        s = char(string(d));
    end
else
    s = '';
end
end

function s = labCodeToString(T, r)
if ismember('LabCode', T.Properties.VariableNames)
    v = T.LabCode(r);
    if iscell(v)
        s = char(string(v{1}));
    elseif isstring(v) || ischar(v)
        s = char(string(v));
    else
        s = char(string(v));
    end
else
    s = '';
end
end
