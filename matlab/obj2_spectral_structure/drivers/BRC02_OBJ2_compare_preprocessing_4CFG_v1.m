%% ============================================================
% BRC02_OBJ2_compare_preprocessing_4CFG_v1
% Objective 2 — Compare 4 preprocessing configurations:
%   (1) RAW
%   (2) SNV
%   (3) SNV + Savitzky–Golay 1st derivative (SG1)
%   (4) SNV + Savitzky–Golay 2nd derivative (SG2)
%
% INPUT
%   - INFILE: Excel table with:
%       * factor columns: Part, Maturity, N2 (names auto-detected)
%       * spectral columns named nm_<integer> (e.g., nm_250 ... nm_1800)
%
% OUTPUT (written under OUTDIR)
%   - Per-configuration folders: OUTDIR/CFG_<cfg>/*
%   - Summary workbook:         OUTDIR/OBJ2_Compare_4CFG.xlsx
%   - Comparison figures:       OUTDIR/COMPARE_FIGS/*
%
% Notes
%   - Spectral plots use Xplot (no autoscaling) to preserve meaningful scale.
%   - PCA / factorial test / ASCA use Xmodel (autoscaled columns) for comparability.
%
% Requirements
%   - MATLAB R2016b+ (local functions in scripts)
%   - Statistics and Machine Learning Toolbox (pca, dummyvar, gscatter, sgolay)
% ============================================================

clear; clc;

% ---- Repository-relative paths (portable) ----
SCRIPT_DIR = fileparts(mfilename('fullpath'));
REPO_ROOT  = fileparts(fileparts(fileparts(SCRIPT_DIR))); % drivers -> obj2 -> matlab -> repo root

%% -------------------------
% USER SETTINGS
%% -------------------------
% Default locations (Zenodo dataset can mirror this structure)
INFILE = fullfile(REPO_ROOT, "data", "processed", "Matriz_Brocoli_SUM_1nm_ASCII.xlsx"); % <-- set your file
OUTDIR = fullfile(REPO_ROOT, "outputs", "obj2", "compare_4cfg");                        % <-- set output dir

if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

% Style
FONT_NAME = 'Times New Roman';
FS_AX  = 10;
FS_LAB = 12;

% PCA export
N_PCS_EXPORT = 10;       % export first N PCs loadings/variance
MAKE_PC13    = true;

% Permutations (heavy). Smoke-test with 199, final with 999+ if needed.
SEED0  = 123;
N_PERM = 499;

% SG parameters
SG_POLY   = 2;
SG_WINDOW = 11;

%% -------------------------
% READ TABLE
%% -------------------------
Tcore = readtable_safely(INFILE);

% Detect factor columns robustly (no accents needed)
colParte = pickVarName(Tcore, ["parte","part"]);
colMad   = pickVarName(Tcore, ["madur","matur"]);
colN2    = pickVarName(Tcore, ["n2","nitro","aplic"]);

assert(colParte ~= "", "Could not detect Part/Parte column.");
assert(colMad   ~= "", "Could not detect Maturity/Maduracion column.");
assert(colN2    ~= "", "Could not detect N2/Nitrogen column.");

%% =========================
% EXTRACT SPECTRA X (nm_250..nm_1800)
% =========================
allNames  = string(Tcore.Properties.VariableNames);
isSpec    = startsWith(allNames, "nm_");
specNames = allNames(isSpec);
assert(~isempty(specNames), 'No spectral columns nm_* found.');

% Parse wavelengths from names and sort
wl = zeros(numel(specNames),1);
for i=1:numel(specNames)
    tok = regexp(specNames(i), '^nm_(\d+)$', 'tokens', 'once');
    assert(~isempty(tok), "Bad spectral column name: %s", specNames(i));
    wl(i) = str2double(tok{1});
end
[wl, ord] = sort(wl, 'ascend');
specNames = specNames(ord);

Xraw = zeros(height(Tcore), numel(specNames));
for j=1:numel(specNames)
    Xraw(:,j) = double(Tcore.(specNames(j)));
end

if any(isnan(Xraw(:)))
    warning('OBJ2:XrawHasNaN', 'Xraw contains NaNs. Check import/export pipeline.');
end

% Check wavelength grid (needed for SG derivatives)
dx = checkUniformGrid(wl);

% Factors as categoricals
Part = categorical(string(Tcore.(colParte)));
Mat  = categorical(string(Tcore.(colMad)));
N2   = categorical(string(Tcore.(colN2)));

%% =========================
% LOOP 4 PREPROCESS CONFIGS
% =========================
cfgList = ["RAW","SNV","SNV_SG1","SNV_SG2"];

% Collectors for comparisons
tblPCA_all    = table();
tblMANOVA_all = table();
tblASCA_all   = table();

% For ASCA effect PC1 loading comparisons (term x cfg)
terms = defineTerms();
nTerms = numel(terms);
effPC1 = cell(nTerms, numel(cfgList));   % each cell: [nWl x 1] loading
effPC1Expl = NaN(nTerms, numel(cfgList));

for c = 1:numel(cfgList)
    cfg = string(cfgList(c));
    fprintf('\n=== OBJ2: Preprocessing config %s ===\n', cfg);

    cfgOut = fullfile(OUTDIR, "CFG_" + cfg);
    if ~exist(cfgOut,'dir'); mkdir(cfgOut); end

    % ---- Apply preprocessing ----
    % Xplot  : preprocessed spectra WITHOUT autoscale (for spectral plots)
    % Xmodel : autoscaled columns (for PCA/MANOVA/ASCA)
    [Xmodel, Xplot, preprocInfo] = preprocessSpectra_4CFG(Xraw, cfg, SG_POLY, SG_WINDOW, dx);

    % ---- Export preprocessing parameters (audit) ----
    save(fullfile(cfgOut, "OBJ2_PreprocInfo_CFG_" + cfg + ".mat"), ...
        "preprocInfo", "wl", "specNames");

    %% ================
    % 2.1 DESCRIPTIVE (on Xplot)
    % ================
    mu = mean(Xplot, 1, 'omitnan');
    sd = std(Xplot, 0, 1, 'omitnan');
    tblMean = table(wl(:), mu(:), sd(:), 'VariableNames', {'Wavelength_nm','Mean','SD'});
    writetable(tblMean, fullfile(cfgOut, "OBJ2_Spectra_MeanSD_CFG_" + cfg + ".xlsx"), 'Sheet', 'MeanSD');

    exportMeanByFactor(Xplot, wl, Part, "Part",     cfgOut, cfg, FONT_NAME, FS_AX, FS_LAB);
    exportMeanByFactor(Xplot, wl, Mat,  "Maturity", cfgOut, cfg, FONT_NAME, FS_AX, FS_LAB);
    exportMeanByFactor(Xplot, wl, N2,   "N2",       cfgOut, cfg, FONT_NAME, FS_AX, FS_LAB);

    % Plot global mean ± SD (Xplot)
    fig = figure('Color','w','Units','pixels','Position',[80 80 1100 520]);
    plot(wl, mu, 'LineWidth', 1.5); hold on;
    plot(wl, mu+sd, '--', 'LineWidth', 1.0);
    plot(wl, mu-sd, '--', 'LineWidth', 1.0);
    xlabel('Wavelength (nm)','FontName',FONT_NAME,'FontSize',FS_LAB);
    ylabel('Preprocessed spectra (a.u.)','FontName',FONT_NAME,'FontSize',FS_LAB);
    title("OBJ2: Global mean \pm SD (CFG " + cfg + ")", 'FontName',FONT_NAME,'FontSize',FS_LAB);
    legend({'Mean','Mean+SD','Mean-SD'}, 'Location','best', 'Box','off');
    set(gca,'FontName',FONT_NAME,'FontSize',FS_AX);
    grid on;
    robustSaveFig(fig, fullfile(cfgOut, "OBJ2_GlobalMeanSD_CFG_" + cfg));

    %% ============
    % 2.2 PCA (on Xmodel)
    % ============
    [coeff, score, ~, ~, explained] = pca(Xmodel, 'Centered', false, 'Algorithm','svd');

    % Export explained variance
    nPC = min(N_PCS_EXPORT, numel(explained));
    tblPCA = table((1:nPC)', explained(1:nPC), cumsum(explained(1:nPC)), ...
        'VariableNames', {'PC','Explained_pct','Cumulative_pct'});
    writetable(tblPCA, fullfile(cfgOut, "OBJ2_PCA_Explained_CFG_" + cfg + ".xlsx"), 'Sheet', 'Explained');

    % Export loadings
    loadTbl = array2table(coeff(:,1:nPC), 'VariableNames', compose("PC%d", 1:nPC));
    loadTbl = addvars(loadTbl, wl(:), 'Before', 1, 'NewVariableNames', 'Wavelength_nm');
    writetable(loadTbl, fullfile(cfgOut, "OBJ2_PCA_Loadings_CFG_" + cfg + ".xlsx"), 'Sheet', 'Loadings');

    % Score plots coloured by each factor separately
    plotPCA_scores(score, explained, Part, "Part",     cfgOut, cfg, FONT_NAME, FS_AX, FS_LAB, MAKE_PC13);
    plotPCA_scores(score, explained, Mat,  "Maturity", cfgOut, cfg, FONT_NAME, FS_AX, FS_LAB, MAKE_PC13);
    plotPCA_scores(score, explained, N2,   "N2",       cfgOut, cfg, FONT_NAME, FS_AX, FS_LAB, MAKE_PC13);

    % Append PCA explained to comparison table
    tmpPCA = tblPCA;
    tmpPCA = addvars(tmpPCA, repmat(cfg,height(tmpPCA),1), 'Before', 1, 'NewVariableNames','CFG');
    tblPCA_all = [tblPCA_all; tmpPCA]; %#ok<AGROW>

    %% =========================================
    % 2.3 MULTIVARIATE FACTORIAL TEST (X space) on Xmodel
    % =========================================
    rng(SEED0);
    results = runFreedmanLaneMANOVA(Xmodel, Part, Mat, N2, terms, N_PERM);

    % BH-FDR across terms within this CFG (6 tests)
    results.q_BH = bh_fdr(results.p_perm);

    writetable(results, fullfile(cfgOut, "OBJ2_SpectralFactorialTest_CFG_" + cfg + ".xlsx"), ...
        'Sheet', 'FactorialTest');

    % Append to comparison table
    tmpM = results;
    tmpM = addvars(tmpM, repmat(cfg,height(tmpM),1), 'Before', 1, 'NewVariableNames','CFG');
    tblMANOVA_all = [tblMANOVA_all; tmpM]; %#ok<AGROW>

    %% ===========================
    % 2.4 ASCA-style decomposition on Xmodel
    % ===========================
    asca = runASCAeffects(Xmodel, Part, Mat, N2, terms);

    writetable(asca.varTable, fullfile(cfgOut, "OBJ2_ASCA_VariancePartition_CFG_" + cfg + ".xlsx"), ...
        'Sheet', 'Variance');

    % Append ASCA variance table
    tmpA = asca.varTable;
    tmpA = addvars(tmpA, repmat(cfg,height(tmpA),1), 'Before', 1, 'NewVariableNames','CFG');
    tblASCA_all = [tblASCA_all; tmpA]; %#ok<AGROW>

    % Effect loadings PC1 per term (and store for cross-CFG comparison)
    for t=1:numel(asca.effects)
        effName = asca.effects(t).name;
        E = asca.effects(t).E;

        [cE, ~, ~, ~, eExpl] = pca(E, 'Centered', false, 'Algorithm','svd');

        effLoad = table(wl(:), cE(:,1), 'VariableNames', {'Wavelength_nm','PC1_loading'});
        writetable(effLoad, fullfile(cfgOut, "OBJ2_ASCA_LoadingPC1_" + effName + "_CFG_" + cfg + ".xlsx"), ...
            'Sheet', 'PC1');

        fig = figure('Color','w','Units','pixels','Position',[80 80 1100 520]);
        plot(wl, cE(:,1), 'LineWidth', 1.4);
        xlabel('Wavelength (nm)','FontName',FONT_NAME,'FontSize',FS_LAB);
        ylabel('PC1 loading (effect PCA)','FontName',FONT_NAME,'FontSize',FS_LAB);
        title("OBJ2 ASCA: " + effName + " effect - PC1 loading (CFG " + cfg + ...
            ", PC1=" + sprintf('%.1f', eExpl(1)) + "%)", ...
            'FontName',FONT_NAME,'FontSize',FS_LAB);
        set(gca,'FontName',FONT_NAME,'FontSize',FS_AX);
        grid on;
        robustSaveFig(fig, fullfile(cfgOut, "OBJ2_ASCA_PC1Loading_" + effName + "_CFG_" + cfg));

        % Store for cross-CFG plot
        termIdx = find(strcmp({terms.name}, effName), 1);
        if ~isempty(termIdx)
            effPC1{termIdx, c} = cE(:,1);
            effPC1Expl(termIdx, c) = eExpl(1);
        end
    end

    fprintf('CFG %s done. Outputs: %s\n', cfg, cfgOut);
end

fprintf('\nOBJ2 complete. Root OUTDIR: %s\n', OUTDIR);

%% =========================
% COMPARISON EXPORTS + FIGS
% =========================
cmpXlsx = fullfile(OUTDIR, "OBJ2_Compare_4CFG.xlsx");
if exist(cmpXlsx,'file'); delete(cmpXlsx); end

writetable(tblPCA_all,    cmpXlsx, 'Sheet','PCA_Explained');
writetable(tblMANOVA_all, cmpXlsx, 'Sheet','FactorialTest');
writetable(tblASCA_all,   cmpXlsx, 'Sheet','ASCA_Variance');

% Wide summaries (easier to compare)
try
    etaWide = unstack(tblMANOVA_all(:,{'CFG','Term','eta2p'}), 'eta2p', 'CFG');
    qWide   = unstack(tblMANOVA_all(:,{'CFG','Term','q_BH'}),   'q_BH',  'CFG');
    writetable(etaWide, cmpXlsx, 'Sheet','eta2p_wide');
    writetable(qWide,   cmpXlsx, 'Sheet','qBH_wide');
catch
    warning('Could not unstack MANOVA tables to wide format. Exported long format only.');
end

% Comparison figures folder
cmpFigDir = fullfile(OUTDIR, "COMPARE_FIGS");
if ~exist(cmpFigDir,'dir'); mkdir(cmpFigDir); end

% (1) PCA explained variance comparison (PC1-3)
fig = figure('Color','w','Units','pixels','Position',[80 80 900 520]);
pcKeep = tblPCA_all.PC <= 3;
Tpc = tblPCA_all(pcKeep,:);
cfgCat = categorical(Tpc.CFG, cfgList, cfgList); %#ok<NASGU>
pcCat  = categorical("PC"+string(Tpc.PC));       %#ok<NASGU>
% Build matrix cfg x pc
M = NaN(numel(cfgList),3);
for i=1:numel(cfgList)
    for p=1:3
        idx = (Tpc.CFG==cfgList(i)) & (Tpc.PC==p);
        if any(idx); M(i,p) = Tpc.Explained_pct(idx); end
    end
end
bar(categorical(cfgList), M);
xlabel('Preprocessing CFG','FontName',FONT_NAME,'FontSize',FS_LAB);
ylabel('Explained variance (%)','FontName',FONT_NAME,'FontSize',FS_LAB);
title('OBJ2 PCA: Explained variance (PC1-3) across CFGs','FontName',FONT_NAME,'FontSize',FS_LAB);
legend({'PC1','PC2','PC3'}, 'Location','best', 'Box','off');
set(gca,'FontName',FONT_NAME,'FontSize',FS_AX);
grid on;
robustSaveFig(fig, fullfile(cmpFigDir, "OBJ2_COMPARE_PCA_Explained_PC1to3"));

% (2) eta2p comparison across CFGs
fig = figure('Color','w','Units','pixels','Position',[80 80 1100 520]);
termOrder = string({terms.name});
etaMat = NaN(numel(termOrder), numel(cfgList));
for t=1:numel(termOrder)
    for c=1:numel(cfgList)
        idx = (tblMANOVA_all.Term==termOrder(t)) & (tblMANOVA_all.CFG==cfgList(c));
        if any(idx); etaMat(t,c) = tblMANOVA_all.eta2p(idx); end
    end
end
bar(categorical(termOrder), etaMat);
xlabel('Term','FontName',FONT_NAME,'FontSize',FS_LAB);
ylabel('Partial eta^2','FontName',FONT_NAME,'FontSize',FS_LAB);
title('OBJ2 Factorial test: partial eta^2 across CFGs','FontName',FONT_NAME,'FontSize',FS_LAB);
legend(cellstr(cfgList), 'Location','best', 'Box','off');
set(gca,'FontName',FONT_NAME,'FontSize',FS_AX);
grid on;
robustSaveFig(fig, fullfile(cmpFigDir, "OBJ2_COMPARE_eta2p_terms"));

% (3) ASCA effect PC1 loading overlays per term
for t=1:nTerms
    % if any empty, skip
    if all(cellfun(@isempty, effPC1(t,:))); continue; end

    fig = figure('Color','w','Units','pixels','Position',[80 80 1100 520]);
    hold on;
    for c=1:numel(cfgList)
        if ~isempty(effPC1{t,c})
            plot(wl, effPC1{t,c}, 'LineWidth', 1.2);
        end
    end
    xlabel('Wavelength (nm)','FontName',FONT_NAME,'FontSize',FS_LAB);
    ylabel('PC1 loading (ASCA effect)','FontName',FONT_NAME,'FontSize',FS_LAB);
    title("OBJ2 ASCA effect PC1 loading overlay: " + termOrder(t), 'FontName',FONT_NAME,'FontSize',FS_LAB);
    legend(cellstr(cfgList), 'Location','best', 'Box','off');
    set(gca,'FontName',FONT_NAME,'FontSize',FS_AX);
    grid on;
    robustSaveFig(fig, fullfile(cmpFigDir, "OBJ2_COMPARE_ASCA_PC1Loading_" + termOrder(t)));
end

fprintf('\nComparison exports done:\n  %s\n  %s\n', cmpXlsx, cmpFigDir);

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function T = readtable_safely(infile)
    % Robust readtable across MATLAB versions
    if ~isfile(infile)
        error('Input file not found: %s', infile);
    end
    try
        T = readtable(infile, 'VariableNamingRule','modify');
    catch
        T = readtable(infile); % fallback
    end
end

function name = pickVarName(T, candidates)
    vars = string(T.Properties.VariableNames);
    name = "";

    % Exact (case-insensitive)
    for c = candidates
        idx = find(strcmpi(vars, c), 1);
        if ~isempty(idx)
            name = vars(idx);
            return;
        end
    end

    % Contains (case-insensitive)
    lvars = lower(vars);
    for c = candidates
        idx = find(contains(lvars, lower(c)), 1);
        if ~isempty(idx)
            name = vars(idx);
            return;
        end
    end
end

function dx = checkUniformGrid(wl)
    d = diff(wl(:));
    dx = median(d);
    if any(abs(d - dx) > 1e-9)
        warning('Wavelength grid is not perfectly uniform. Using median dx=%.6g for SG derivatives.', dx);
    end
end

function [Xmodel, Xplot, info] = preprocessSpectra_4CFG(Xraw, cfg, sgPoly, sgWindow, dx)
    cfg = upper(string(cfg));
    info = struct();
    info.cfg = cfg;
    info.sg_poly   = sgPoly;
    info.sg_window = sgWindow;
    info.dx        = dx;

    switch cfg
        case "RAW"
            Xwork = Xraw;
            info.SNV = false;
            info.sg_derivative_order = 0;

        case "SNV"
            Xwork = snv_rows(Xraw);
            info.SNV = true;
            info.sg_derivative_order = 0;

        case "SNV_SG1"
            Xsnv  = snv_rows(Xraw);
            Xwork = sg_derivative(Xsnv, 1, sgPoly, sgWindow, dx);
            info.SNV = true;
            info.sg_derivative_order = 1;

        case "SNV_SG2"
            Xsnv  = snv_rows(Xraw);
            Xwork = sg_derivative(Xsnv, 2, sgPoly, sgWindow, dx);
            info.SNV = true;
            info.sg_derivative_order = 2;

        otherwise
            error('Unknown cfg: %s', cfg);
    end

    % For spectral plots (meaningful scale)
    Xplot = Xwork;

    % For PCA/MANOVA/ASCA (comparable feature scaling)
    [Xmodel, mu, sigma] = autoscale_cols(Xwork);
    info.colMean = mu;
    info.colSD   = sigma;
end

function Xs = snv_rows(X)
    mu = mean(X, 2, 'omitnan');
    sd = std(X, 0, 2, 'omitnan');
    sd(sd==0) = 1;
    Xs = (X - mu) ./ sd;
end

function [Xa, mu, sd] = autoscale_cols(X)
    mu = mean(X, 1, 'omitnan');
    sd = std(X, 0, 1, 'omitnan');
    sd(sd==0) = 1;
    Xa = (X - mu) ./ sd;
end

function Xd = sg_derivative(X, derivOrder, polyOrder, frameLen, dx)
    if mod(frameLen,2) == 0
        error('Savitzky-Golay frameLen must be odd.');
    end
    if derivOrder ~= 1 && derivOrder ~= 2
        error('Only 1st or 2nd derivative supported (derivOrder=1 or 2).');
    end

    [~, G] = sgolay(polyOrder, frameLen);

    % Derivative filter coefficients
    % order d => column (d+1)
    h = factorial(derivOrder) / (dx^derivOrder) * G(:, derivOrder+1);

    half = (frameLen-1)/2;
    Xd = zeros(size(X));

    for i=1:size(X,1)
        x = X(i,:);
        % mirror padding
        left  = x(half:-1:1);
        right = x(end:-1:end-half+1);
        xpad = [left, x, right];

        y = conv(xpad, flipud(h), 'valid'); % length == size(x,2)
        Xd(i,:) = y;
    end
end

function exportMeanByFactor(Xplot, wl, g, gName, outDir, cfg, fontName, fsAx, fsLab)
    lv = categories(g);
    M = NaN(numel(lv), numel(wl));
    for i=1:numel(lv)
        idx = (g == lv{i});
        M(i,:) = mean(Xplot(idx,:), 1, 'omitnan');
    end

    tbl = array2table(M, 'VariableNames', compose("nm_%d", wl));
    tbl = addvars(tbl, string(lv(:)), 'Before', 1, 'NewVariableNames', gName);
    writetable(tbl, fullfile(outDir, "OBJ2_MeanSpectraBy_" + gName + "_CFG_" + cfg + ".xlsx"), ...
        'Sheet', 'Means');

    fig = figure('Color','w','Units','pixels','Position',[80 80 1100 520]);
    hold on;
    for i=1:numel(lv)
        plot(wl, M(i,:), 'LineWidth', 1.2);
    end
    xlabel('Wavelength (nm)','FontName',fontName,'FontSize',fsLab);
    ylabel('Preprocessed spectra (a.u.)','FontName',fontName,'FontSize',fsLab);
    title("OBJ2: Mean spectra by " + gName + " (CFG " + cfg + ")", ...
        'FontName',fontName,'FontSize',fsLab);
    legend(string(lv), 'Location','best', 'Box','off');
    set(gca,'FontName',fontName,'FontSize',fsAx);
    grid on;
    robustSaveFig(fig, fullfile(outDir, "OBJ2_MeanSpectraBy_" + gName + "_CFG_" + cfg));
end

function plotPCA_scores(score, explained, g, gName, outDir, cfg, fontName, fsAx, fsLab, makePC13)
    fig = figure('Color','w','Units','pixels','Position',[80 80 820 650]);
    gscatter(score(:,1), score(:,2), g);
    xlabel(sprintf('PC1 (%.1f%%)', explained(1)), 'FontName',fontName,'FontSize',fsLab);
    ylabel(sprintf('PC2 (%.1f%%)', explained(2)), 'FontName',fontName,'FontSize',fsLab);
    title("OBJ2 PCA scores by " + gName + " (CFG " + cfg + ")", 'FontName',fontName,'FontSize',fsLab);
    legend('Location','best', 'Box','off');
    set(gca,'FontName',fontName,'FontSize',fsAx);
    grid on;
    robustSaveFig(fig, fullfile(outDir, "OBJ2_PCA_PC1PC2_by_" + gName + "_CFG_" + cfg));

    if makePC13
        fig = figure('Color','w','Units','pixels','Position',[80 80 820 650]);
        gscatter(score(:,1), score(:,3), g);
        xlabel(sprintf('PC1 (%.1f%%)', explained(1)), 'FontName',fontName,'FontSize',fsLab);
        ylabel(sprintf('PC3 (%.1f%%)', explained(3)), 'FontName',fontName,'FontSize',fsLab);
        title("OBJ2 PCA scores by " + gName + " (PC1-PC3, CFG " + cfg + ")", 'FontName',fontName,'FontSize',fsLab);
        legend('Location','best', 'Box','off');
        set(gca,'FontName',fontName,'FontSize',fsAx);
        grid on;
        robustSaveFig(fig, fullfile(outDir, "OBJ2_PCA_PC1PC3_by_" + gName + "_CFG_" + cfg));
    end
end

function robustSaveFig(fig, basePathNoExt)
    drawnow;
    try
        saveas(fig, basePathNoExt + ".fig");
    catch
        warning('SaveFIG failed for: %s', basePathNoExt);
    end
    try
        exportgraphics(fig, basePathNoExt + ".png", 'Resolution', 400);
    catch
        print(fig, basePathNoExt + ".png", '-dpng', '-r400');
    end
    close(fig);
end

function terms = defineTerms()
    terms = struct();
    terms(1).name = "Part";
    terms(1).type = "main";
    terms(2).name = "Maturity";
    terms(2).type = "main";
    terms(3).name = "N2";
    terms(3).type = "main";
    % Use ASCII-safe interaction labels (portable filenames / encoding)
    terms(4).name = "PartxMaturity";
    terms(4).type = "int";
    terms(5).name = "PartxN2";
    terms(5).type = "int";
    terms(6).name = "MaturityxN2";
    terms(6).type = "int";
end

function Z = buildDesign(Part, Mat, N2, include)
    n = numel(Part);
    Z = ones(n,1);

    Dp = dummyvar(Part); Dp = Dp(:,2:end);
    Dm = dummyvar(Mat);  Dm = Dm(:,2:end);
    Dn = dummyvar(N2);   Dn = Dn(:,2:end);

    if include.includePart; Z = [Z, Dp]; end %#ok<AGROW>
    if include.includeMat;  Z = [Z, Dm]; end %#ok<AGROW>
    if include.includeN2;   Z = [Z, Dn]; end %#ok<AGROW>

    if include.includePM; Z = [Z, interactionDummy(Dp, Dm)]; end %#ok<AGROW>
    if include.includePN; Z = [Z, interactionDummy(Dp, Dn)]; end %#ok<AGROW>
    if include.includeMN; Z = [Z, interactionDummy(Dm, Dn)]; end %#ok<AGROW>
end

function Dint = interactionDummy(D1, D2)
    [n, p1] = size(D1);
    p2 = size(D2,2);
    Dint = zeros(n, p1*p2);
    k = 1;
    for i=1:p1
        for j=1:p2
            Dint(:,k) = D1(:,i) .* D2(:,j);
            k = k + 1;
        end
    end
end

function resTbl = runFreedmanLaneMANOVA(X, Part, Mat, N2, terms, nPerm)
    n = size(X,1);

    incFull = struct('includePart',true,'includeMat',true,'includeN2',true, ...
                     'includePM',true,'includePN',true,'includeMN',true);

    Zfull = buildDesign(Part, Mat, N2, incFull);
    dfFull = rank(Zfull);
    dfErr  = n - dfFull;

    XhatFull = Zfull * (Zfull \ X);
    SSEfull  = sum((X - XhatFull).^2, 'all');
    SST      = sum(X.^2, 'all');

    resTbl = table('Size',[numel(terms) 10], ...
        'VariableTypes', {'string','double','double','double','double','double','double','double','double','double'}, ...
        'VariableNames', {'Term','df_term','df_error','SS_term','SS_error','F','p_perm','R2','eta2p','nPerm'});

    for t=1:numel(terms)
        termName = terms(t).name;

        incRed = incFull;
        switch termName
            case "Part";            incRed.includePart = false;
            case "Maturity";        incRed.includeMat  = false;
            case "N2";              incRed.includeN2   = false;
            case "PartxMaturity";   incRed.includePM   = false;
            case "PartxN2";         incRed.includePN   = false;
            case "MaturityxN2";     incRed.includeMN   = false;
            otherwise
                error('Unknown term: %s', termName);
        end

        Zred = buildDesign(Part, Mat, N2, incRed);
        dfRed  = rank(Zred);
        dfTerm = dfFull - dfRed;

        XhatRed = Zred * (Zred \ X);
        SSterm_obs = sum((XhatFull - XhatRed).^2, 'all');

        Fobs = (SSterm_obs/dfTerm) / (SSEfull/dfErr);

        % Freedman-Lane permutations
        XhatRed0 = XhatRed;
        Ered0    = X - XhatRed0;

        Fperm = zeros(nPerm,1);
        for b=1:nPerm
            idx = randperm(n);
            Xb  = XhatRed0 + Ered0(idx,:);

            XhatFull_b = Zfull * (Zfull \ Xb);
            SSEfull_b  = sum((Xb - XhatFull_b).^2, 'all');

            XhatRed_b  = Zred * (Zred \ Xb);
            SSterm_b   = sum((XhatFull_b - XhatRed_b).^2, 'all');

            Fperm(b) = (SSterm_b/dfTerm) / (SSEfull_b/dfErr);
        end

        pPerm = (1 + sum(Fperm >= Fobs)) / (1 + nPerm);

        resTbl.Term(t)     = termName;
        resTbl.df_term(t)  = dfTerm;
        resTbl.df_error(t) = dfErr;
        resTbl.SS_term(t)  = SSterm_obs;
        resTbl.SS_error(t) = SSEfull;
        resTbl.F(t)        = Fobs;
        resTbl.p_perm(t)   = pPerm;
        resTbl.R2(t)       = SSterm_obs / SST;
        resTbl.eta2p(t)    = SSterm_obs / (SSterm_obs + SSEfull);
        resTbl.nPerm(t)    = nPerm;

        fprintf('Term %-14s | F=%.3f | p_perm=%.4g | R2=%.3f\n', termName, Fobs, pPerm, resTbl.R2(t));
    end
end

function out = runASCAeffects(X, Part, Mat, N2, terms)
    incFull = struct('includePart',true,'includeMat',true,'includeN2',true, ...
                     'includePM',true,'includePN',true,'includeMN',true);

    Zfull = buildDesign(Part, Mat, N2, incFull);
    XhatFull = Zfull * (Zfull \ X);
    Efull = X - XhatFull;

    SST = sum(X.^2, 'all');
    SSE = sum(Efull.^2, 'all');

    varTable = table('Size',[numel(terms)+1 4], ...
        'VariableTypes', {'string','double','double','double'}, ...
        'VariableNames', {'Component','SS','Pct_of_Total','Pct_of_ModelPlusError'});

    effects = struct('name',{},'E',{});

    ssSum = 0;
    for t=1:numel(terms)
        termName = terms(t).name;

        incRed = incFull;
        switch termName
            case "Part";            incRed.includePart = false;
            case "Maturity";        incRed.includeMat  = false;
            case "N2";              incRed.includeN2   = false;
            case "PartxMaturity";   incRed.includePM   = false;
            case "PartxN2";         incRed.includePN   = false;
            case "MaturityxN2";     incRed.includeMN   = false;
        end

        Zred = buildDesign(Part, Mat, N2, incRed);
        XhatRed = Zred * (Zred \ X);

        Eff = XhatFull - XhatRed;
        SS_eff = sum(Eff.^2, 'all');
        ssSum = ssSum + SS_eff;

        effects(end+1).name = termName; %#ok<AGROW>
        effects(end).E = Eff;

        varTable.Component(t) = termName;
        varTable.SS(t) = SS_eff;
        varTable.Pct_of_Total(t) = 100 * SS_eff / SST;
        varTable.Pct_of_ModelPlusError(t) = 100 * SS_eff / (ssSum + SSE);
    end

    varTable.Component(end) = "Error";
    varTable.SS(end) = SSE;
    varTable.Pct_of_Total(end) = 100 * SSE / SST;
    varTable.Pct_of_ModelPlusError(end) = 100 * SSE / (ssSum + SSE);

    out = struct();
    out.varTable = varTable;
    out.effects  = effects;
end

function q = bh_fdr(p)
    p = p(:);
    q = NaN(size(p));
    ok = ~isnan(p);
    p0 = p(ok);

    [ps, idx] = sort(p0, 'ascend');
    m = numel(ps);
    if m==0; return; end

    qs = ps .* (m ./ (1:m)');
    for i=m-1:-1:1
        qs(i) = min(qs(i), qs(i+1));
    end
    qs(qs>1) = 1;

    q0 = NaN(size(p0));
    q0(idx) = qs;
    q(ok) = q0;
end
