% BRC02_Fig23_SpectraPCA_SNV_SG1st_SG2nd_v1.m
% -------------------------------------------------------------------------
% Manuscript figures for Section 3.2 (main text): two double-panel figures.
%
% Fig. 2 — Mean spectra by plant part (Part)
%   (A) SNV
%   (B) SNV + Savitzky–Golay 1st derivative
%
% Fig. 3 — PCA scores (PC1 vs PC2) coloured by plant part (Part)
%   (A) SNV
%   (B) SNV + Savitzky–Golay 2nd derivative
%
% Input:
%   - Matriz_Brocoli_SUM_1nm_ASCII.xlsx
%   - Spectral columns must be named nm_<integer> (e.g., nm_250 ... nm_1800)
%
% Method notes (consistent with OBJ2 compare-preprocessing logic):
%   - Mean spectra are plotted by Part only (one curve per Part level).
%   - PCA is performed on autoscaled X (column-wise mean-centre + unit-variance).
%   - Savitzky–Golay parameters are user-defined (default: poly=2, window=11, delta=1 nm).
% -------------------------------------------------------------------------

clear; clc;

%% =========================
% USER PARAMETERS (EDIT)
% =========================
INPUT_XLSX  = 'Matriz_Brocoli_SUM_1nm_ASCII.xlsx';
SHEET_NAME  = 'Matriz';

OUTDIR = fullfile(pwd, 'FIG_3p2_OUT');
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

% Optional subset filter (use only if Section 3.2 is intended for a specific core)
USE_CORE_FILTER = true;
CULTIVAR_CORE   = "Pathernon";    % set to the exact label used in the spreadsheet
EXTRACTION_CORE = "Ultrasonido";  % set to the exact label used in the spreadsheet

% Savitzky–Golay parameters
SG_POLY_ORDER = 2;
SG_FRAME_LEN  = 11;   % must be odd
DELTA_NM      = 1;

% Plot style
FONT_NAME = 'Times New Roman';
FS_AX  = 10;
FS_LAB = 12;

%% =========================
% LOAD DATA (PRESERVE HEADERS)
% =========================
fprintf('Loading: %s\n', INPUT_XLSX);
opts = detectImportOptions(INPUT_XLSX, 'Sheet', SHEET_NAME);
opts.VariableNamingRule = 'preserve';
T = readtable(INPUT_XLSX, opts);

fprintf('Rows: %d | Cols: %d\n', height(T), width(T));

%% =========================
% IDENTIFY REQUIRED COLUMNS
% =========================
colVariedad = pickVarName(T, ["Variedad","Variety","Cultivar"]);
colParte    = pickVarName(T, ["Parte","Part"]);
colExt      = pickVarName(T, ["Extraccion","Extracción","Extraction"]);

assert(colParte~="", 'Missing Part/Parte column (check matrix headers).');

if USE_CORE_FILTER
    assert(colVariedad~="" && colExt~="", ...
        'Missing Cultivar/Variety and/or Extraction column required for core filtering.');
end

%% =========================
% OPTIONAL FILTER (CORE SUBSET)
% =========================
if USE_CORE_FILTER
    varietyStr = string(T.(colVariedad));
    extrStr    = string(T.(colExt));
    isCore = strcmpi(strtrim(varietyStr), CULTIVAR_CORE) & strcmpi(strtrim(extrStr), EXTRACTION_CORE);
    T = T(isCore, :);
    fprintf('Core subset (Cultivar=%s, Extraction=%s): n = %d\n', CULTIVAR_CORE, EXTRACTION_CORE, height(T));
else
    fprintf('No core filter applied: n = %d\n', height(T));
end

assert(height(T) > 5, 'Too few rows after filtering. Check core filter settings.');

%% =========================
% EXTRACT SPECTRA X (nm_*)
% =========================
allNames  = string(T.Properties.VariableNames);
isSpec    = startsWith(allNames, "nm_");
specNames = allNames(isSpec);
assert(~isempty(specNames), 'No spectral columns found (expected "nm_*").');

% Parse wavelengths and sort
wl = nan(numel(specNames),1);
for i = 1:numel(specNames)
    tok = regexp(specNames(i), '^nm_(\d+)$', 'tokens', 'once');
    if ~isempty(tok); wl(i) = str2double(tok{1}); end
end
[wl, ord] = sort(wl, 'ascend');
specNames = specNames(ord);

Xraw = double(table2array(T(:, specNames)));

if any(isnan(Xraw(:)))
    warning('Xraw contains NaNs. Check import/export pipeline.');
end

Part = categorical(string(T.(colParte)));
Part = removecats(Part);

%% =========================
% PREPROCESSING: SNV, SNV+SG 1st, SNV+SG 2nd
% =========================
X_snv      = snv_rows(Xraw);
X_snv_sg1  = sg_derivative_rows(X_snv, SG_POLY_ORDER, SG_FRAME_LEN, 1, DELTA_NM);
X_snv_sg2  = sg_derivative_rows(X_snv, SG_POLY_ORDER, SG_FRAME_LEN, 2, DELTA_NM);

%% =========================
% FIGURE 2: MEAN SPECTRA BY PART (SNV vs SG 1st)
% =========================
fig = figure('Color','w','Units','pixels','Position',[80 80 1200 520]);
tlo = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

% (A) SNV
nexttile; hold on;
plotMeanByGroup(wl, X_snv, Part, 1.4);
xlabel('Wavelength (nm)', 'FontName',FONT_NAME, 'FontSize',FS_LAB);
ylabel('SNV-normalised reflectance (a.u.)', 'FontName',FONT_NAME, 'FontSize',FS_LAB);
title('(A) SNV', 'FontName',FONT_NAME, 'FontSize',FS_LAB);
set(gca,'FontName',FONT_NAME,'FontSize',FS_AX); grid on;
legend(categories(Part), 'Location','best', 'Box','off');

% (B) SNV + SG 1st derivative
nexttile; hold on;
plotMeanByGroup(wl, X_snv_sg1, Part, 1.4);
xlabel('Wavelength (nm)', 'FontName',FONT_NAME, 'FontSize',FS_LAB);
ylabel('SNV + SG 1st derivative (a.u.)', 'FontName',FONT_NAME, 'FontSize',FS_LAB);
title('(B) SNV + SG 1st derivative', 'FontName',FONT_NAME, 'FontSize',FS_LAB);
set(gca,'FontName',FONT_NAME,'FontSize',FS_AX); grid on;
legend(categories(Part), 'Location','best', 'Box','off');

title(tlo, 'Mean spectra by plant part (visual reference)', 'FontName',FONT_NAME, 'FontSize',FS_LAB);

robustSaveFig(fig, fullfile(OUTDIR, 'Fig2_Spectra_SNV_vs_SG1st'));

%% =========================
% FIGURE 3: PCA SCORES BY PART (SNV vs SG 2nd)
% =========================
% Autoscale (column-wise) before PCA
[X_snv_z, ~, ~]     = autoscale_cols(X_snv);
[X_snv_sg2_z, ~, ~] = autoscale_cols(X_snv_sg2);

[~, score_snv, ~, ~, expl_snv] = pca(X_snv_z, 'Centered', false, 'Algorithm','svd');
[~, score_sg2, ~, ~, expl_sg2] = pca(X_snv_sg2_z, 'Centered', false, 'Algorithm','svd');

fig = figure('Color','w','Units','pixels','Position',[80 80 1200 520]);
tlo = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');

% (A) PCA under SNV
nexttile;
gscatter(score_snv(:,1), score_snv(:,2), Part);
xlabel(sprintf('PC1 (%.1f%%)', expl_snv(1)), 'FontName',FONT_NAME, 'FontSize',FS_LAB);
ylabel(sprintf('PC2 (%.1f%%)', expl_snv(2)), 'FontName',FONT_NAME, 'FontSize',FS_LAB);
title('(A) PCA scores (SNV)', 'FontName',FONT_NAME, 'FontSize',FS_LAB);
set(gca,'FontName',FONT_NAME,'FontSize',FS_AX); grid on;
legend('Location','best', 'Box','off');

% (B) PCA under SNV + SG 2nd derivative
nexttile;
gscatter(score_sg2(:,1), score_sg2(:,2), Part);
xlabel(sprintf('PC1 (%.1f%%)', expl_sg2(1)), 'FontName',FONT_NAME, 'FontSize',FS_LAB);
ylabel(sprintf('PC2 (%.1f%%)', expl_sg2(2)), 'FontName',FONT_NAME, 'FontSize',FS_LAB);
title('(B) PCA scores (SNV + SG 2nd derivative)', 'FontName',FONT_NAME, 'FontSize',FS_LAB);
set(gca,'FontName',FONT_NAME,'FontSize',FS_AX); grid on;
legend('Location','best', 'Box','off');

title(tlo, 'PCA score space by plant part', 'FontName',FONT_NAME, 'FontSize',FS_LAB);

robustSaveFig(fig, fullfile(OUTDIR, 'Fig3_PCA_SNV_vs_SG2nd'));

fprintf('\nDone.\nOutputs in: %s\n', OUTDIR);

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function name = pickVarName(T, candidates)
    vars = string(T.Properties.VariableNames);
    name = "";
    % exact match (case-insensitive)
    for c = candidates
        idx = find(strcmpi(vars, c), 1);
        if ~isempty(idx); name = vars(idx); return; end
    end
    % partial match (case-insensitive)
    for c = candidates
        idx = find(contains(lower(vars), lower(c)), 1);
        if ~isempty(idx); name = vars(idx); return; end
    end
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

function Xd = sg_derivative_rows(X, polyOrder, frameLen, derivOrder, delta)
    if mod(frameLen,2) == 0
        error('Savitzky-Golay frameLen must be odd.');
    end
    [~, G] = sgolay(polyOrder, frameLen);

    % Derivative filter coefficients (column derivOrder+1 in G)
    h = factorial(derivOrder) / (delta^derivOrder) * G(:, derivOrder+1);

    half = (frameLen-1)/2;
    Xd = zeros(size(X));

    for i = 1:size(X,1)
        x = X(i,:);

        % mirror padding (reduce endpoint artefacts)
        left  = x(half:-1:1);
        right = x(end:-1:end-half+1);
        xpad  = [left, x, right];

        y = conv(xpad, flipud(h), 'valid'); % length == numel(x)
        Xd(i,:) = y;
    end
end

function plotMeanByGroup(wl, X, g, lw)
    lv = categories(g);
    for i = 1:numel(lv)
        idx = (g == lv{i});
        mu = mean(X(idx,:), 1, 'omitnan');
        plot(wl, mu, 'LineWidth', lw);
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
