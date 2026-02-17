function O2_SpectralExploration_obj2_run_v3(inputXlsx, sheetName)
% O2_SpectralExploration_obj2_run_v3
% -------------------------------------------------------------------------
% Objective 2 — Hyperspectral exploration (PCA) — Q1-ready outputs
%
% DELIVERABLES (audit-ready)
%   Figures (PNG 300 dpi + MATLAB FIG):
%     Main (aligned to manuscript Fig. 3A–D):
%       - Fig3A_FX10_raw_byDate
%       - Fig3B_FX10_SNV_byIrrigation
%       - Fig3C_FX17_raw_byDate
%       - Fig3D_FX17_SNV_byIrrigation
%
%     Loadings (aligned to manuscript Fig. 4 A–D):
%       - Loadings_FX10_RAW
%       - Loadings_FX10_SNV
%       - Loadings_FX17_RAW
%       - Loadings_FX17_SNV
%
%     Mean + difference spectra by irrigation (aligned to manuscript Fig. 5 A–D):
%       - MeanDiffSpectra_FX10_RAW_byIrrigation
%       - MeanDiffSpectra_FX10_SNV_byIrrigation
%       - MeanDiffSpectra_FX17_RAW_byIrrigation
%       - MeanDiffSpectra_FX17_SNV_byIrrigation
%
%     Supplementary diagnostics (optional; recommended for audit trail):
%       - FigS_FX10_RAW_byIrrigation
%       - FigS_FX10_SNV_byDate
%       - FigS_FX17_RAW_byIrrigation
%       - FigS_FX17_SNV_byDate
%
%   Excel summary:
%     - Explained (PC1–PC5 variance %)
%     - ANOVA_PC1_PC2 (p-values + eta^2 for SamplingDate and IrrigationRegime)
%     - TopLoadings (Top-10 |loading| wavelengths for PC1 and PC2)
%
% OUTPUT LOCATION (strict):
%   All outputs are written into a folder named "Objetivo_2" placed in the
%   root directory of the input spreadsheet.
%
% USAGE:
%   O2_SpectralExploration_obj2_run_v3('Matriz_CHEM_HSI_MASTER_96.xlsx','Matriz');
%
% NOTES:
%   - PCA is computed via SVD on mean-centred data (p >> n robust).
%   - SNV is applied row-wise (per sample).
%   - No Hotelling T^2 is computed (avoids collinearity warnings; equivalent PCA).
% -------------------------------------------------------------------------

    if nargin < 1 || isempty(inputXlsx)
        inputXlsx = 'Matriz_CHEM_HSI_MASTER_96.xlsx';
    end
    if nargin < 2 || isempty(sheetName)
        sheetName = 'Matriz';
    end

    % ---------- Output folder: <root>/Objetivo_2 ----------
    rootDir = fileparts(inputXlsx);
    if isempty(rootDir)
        rootDir = pwd;
    end
    outDir = fullfile(rootDir, 'Objetivo_2');
    if ~exist(outDir,'dir'); mkdir(outDir); end

    % Run tag for log/Excel filenames (prevents overwriting; still inside Objetivo_2)
    runTag = datestr(now,'yyyymmdd_HHMMSS');

    logFile = fullfile(outDir, ['O2_log_' runTag '.txt']);
    diary(logFile); diary on;
    fprintf('=== Objective 2 | Spectral exploration (PCA) | v3 ===\n');
    fprintf('Input: %s | Sheet: %s\n', inputXlsx, sheetName);
    fprintf('Output folder: %s\n', outDir);
    fprintf('Run tag: %s\n', runTag);

    % ---------- Force a white theme (robust against MATLAB defaults) ----------
    set(groot,'defaultFigureColor','w');
    set(groot,'defaultAxesColor','w');
    set(groot,'defaultAxesXColor','k');
    set(groot,'defaultAxesYColor','k');
    set(groot,'defaultTextColor','k');

    % ---------- Load master matrix ----------
    opts = detectImportOptions(inputXlsx, 'Sheet', sheetName);
    opts.VariableNamingRule = 'preserve';
    T = readtable(inputXlsx, opts);

    % ---------- Required metadata columns ----------
    required = {'SamplingDate','IrrigationRegime','VineyardFloorManagement'};
    for k = 1:numel(required)
        if ~ismember(required{k}, T.Properties.VariableNames)
            error('Missing required column: %s (check the matrix).', required{k});
        end
    end

    % ---------- Clean / enforce metadata categoricals ----------
    [DateCat, IrrCat, MngCat] = enforceCategoricals_obj2(T); %#ok<NASGU>

    % ---------- Extract spectra blocks ----------
    [X10, wl10] = extractSpectraPrefix(T, 'FX10_nm_');  % VIS–NIR
    [X17, wl17] = extractSpectraPrefix(T, 'FX17_nm_');  % SWIR

    % ---------- Settings ----------
    NCOMP = 5;
    DO_LOADINGS = true;
    DO_MEAN_DIFF = true;
    DO_SUPP_DIAGNOSTICS = true;

    % ---------- PCA jobs ----------
    jobs = {
        'FX10','RAW',X10,wl10,DateCat,IrrCat;
        'FX10','SNV',X10,wl10,DateCat,IrrCat;
        'FX17','RAW',X17,wl17,DateCat,IrrCat;
        'FX17','SNV',X17,wl17,DateCat,IrrCat
    };

    Summary = {};
    SummaryHeader = {'Sensor','Preproc','nSamples','nWavelengths','PC1_%','PC2_%','PC3_%','PC4_%','PC5_%'};

    ANOVA = {};
    ANOVAHeader = {'Sensor','Preproc','Factor','p_PC1','eta2_PC1','p_PC2','eta2_PC2','nGroups'};

    TopLoad = {};
    TopLoadHeader = {'Sensor','Preproc','PC','Rank','Wavelength_nm','Loading'};

    for i = 1:size(jobs,1)
        sensor = jobs{i,1};
        pp     = jobs{i,2};
        X      = jobs{i,3};
        wl     = jobs{i,4};
        DateC  = jobs{i,5};
        IrrC   = jobs{i,6};

        if isempty(X)
            fprintf('[WARN] %s %s: empty spectra matrix, skipped.\n', sensor, pp);
            continue;
        end

        fprintf('\n--- %s (%s) ---\n', sensor, pp);

        % QC: remove rows with any NaN/Inf in spectra
        badRow = any(~isfinite(X),2);
        if any(badRow)
            fprintf('[QC] Removing %d rows with NaN/Inf in spectra.\n', sum(badRow));
            X = X(~badRow,:);
            DateC = DateC(~badRow);
            IrrC  = IrrC(~badRow);
        end

        % QC: remove near-constant wavelengths
        v = var(X,0,1);
        keep = v > 1e-12;
        if any(~keep)
            fprintf('[QC] Removing %d near-constant wavelengths.\n', sum(~keep));
            X = X(:,keep);
            wl = wl(keep);
        end

        % Preprocess
        if strcmpi(pp,'SNV')
            Xp = snvRows(X);
        else
            Xp = X;
        end

        % PCA via SVD (mean-centred)
        [coeff, score, explained] = pcaSVD(Xp, NCOMP);

        % Summary row
        pc = nan(1,5);
        pc(1:numel(explained)) = explained(:)';
        Summary(end+1,:) = {sensor,pp,size(Xp,1),size(Xp,2),pc(1),pc(2),pc(3),pc(4),pc(5)}; %#ok<AGROW>

        % --- Main score plots (Fig. 3A–D) ---
        if strcmp(sensor,'FX10') && strcmpi(pp,'RAW')
            f = plotScores_obj2(score, explained, DateC, 'Sampling date', sprintf('%s raw reflectance',sensor));
            exportPub_obj2(f, outDir, 'Fig3A_FX10_raw_byDate');
        elseif strcmp(sensor,'FX10') && strcmpi(pp,'SNV')
            f = plotScores_obj2(score, explained, IrrC, 'Irrigation regime', sprintf('%s SNV',sensor));
            exportPub_obj2(f, outDir, 'Fig3B_FX10_SNV_byIrrigation');
        elseif strcmp(sensor,'FX17') && strcmpi(pp,'RAW')
            f = plotScores_obj2(score, explained, DateC, 'Sampling date', sprintf('%s raw reflectance',sensor));
            exportPub_obj2(f, outDir, 'Fig3C_FX17_raw_byDate');
        elseif strcmp(sensor,'FX17') && strcmpi(pp,'SNV')
            f = plotScores_obj2(score, explained, IrrC, 'Irrigation regime', sprintf('%s SNV',sensor));
            exportPub_obj2(f, outDir, 'Fig3D_FX17_SNV_byIrrigation');
        end

        % --- Supplementary diagnostics (recommended) ---
        if DO_SUPP_DIAGNOSTICS
            if strcmp(sensor,'FX10') && strcmpi(pp,'RAW')
                fS = plotScores_obj2(score, explained, IrrC, 'Irrigation regime', sprintf('%s raw reflectance',sensor));
                exportPub_obj2(fS, outDir, 'FigS_FX10_RAW_byIrrigation');
            elseif strcmp(sensor,'FX10') && strcmpi(pp,'SNV')
                fS = plotScores_obj2(score, explained, DateC, 'Sampling date', sprintf('%s SNV',sensor));
                exportPub_obj2(fS, outDir, 'FigS_FX10_SNV_byDate');
            elseif strcmp(sensor,'FX17') && strcmpi(pp,'RAW')
                fS = plotScores_obj2(score, explained, IrrC, 'Irrigation regime', sprintf('%s raw reflectance',sensor));
                exportPub_obj2(fS, outDir, 'FigS_FX17_RAW_byIrrigation');
            elseif strcmp(sensor,'FX17') && strcmpi(pp,'SNV')
                fS = plotScores_obj2(score, explained, DateC, 'Sampling date', sprintf('%s SNV',sensor));
                exportPub_obj2(fS, outDir, 'FigS_FX17_SNV_byDate');
            end
        end

        % --- Optional: loadings (Fig. 4) ---
        if DO_LOADINGS
            fL = plotLoadings_obj2(wl, coeff, explained, sprintf('%s %s', sensor, pp));
            exportPub_obj2(fL, outDir, sprintf('Loadings_%s_%s', sensor, pp));
        end

        % --- Optional: mean + difference spectra by irrigation (Fig. 5) ---
        if DO_MEAN_DIFF
            fM = plotMeanAndDiffSpectra_obj2(wl, Xp, IrrC, sprintf('%s %s', sensor, pp));
            exportPub_obj2(fM, outDir, sprintf('MeanDiffSpectra_%s_%s_byIrrigation', sensor, pp));
        end

        % --- ANOVA on PC1/PC2 for SamplingDate and IrrigationRegime (diagnostic) ---
        [p1,eta1,ng] = anova1_p_eta(score(:,1), DateC);
        [p2,eta2,~ ] = anova1_p_eta(score(:,2), DateC);
        ANOVA(end+1,:) = {sensor,pp,'SamplingDate',p1,eta1,p2,eta2,ng}; %#ok<AGROW>

        [p1,eta1,ng] = anova1_p_eta(score(:,1), IrrC);
        [p2,eta2,~ ] = anova1_p_eta(score(:,2), IrrC);
        ANOVA(end+1,:) = {sensor,pp,'IrrigationRegime',p1,eta1,p2,eta2,ng}; %#ok<AGROW>

        % --- Top loadings wavelengths (PC1 & PC2) ---
        for pcIdx = 1:2
            [~,ord] = sort(abs(coeff(:,pcIdx)),'descend');
            topK = min(10,numel(ord));
            for r = 1:topK
                j = ord(r);
                TopLoad(end+1,:) = {sensor,pp,sprintf('PC%d',pcIdx),r,wl(j),coeff(j,pcIdx)}; %#ok<AGROW>
            end
        end
    end

    % ---------- Write Excel summary (inside Objetivo_2) ----------
    outXlsx = fullfile(outDir, ['O2_PCA_Summary_' runTag '.xlsx']);
    writecell(SummaryHeader, outXlsx, 'Sheet','Explained', 'Range','A1');
    writecell(Summary,       outXlsx, 'Sheet','Explained', 'Range','A2');

    writecell(ANOVAHeader, outXlsx, 'Sheet','ANOVA_PC1_PC2', 'Range','A1');
    writecell(ANOVA,       outXlsx, 'Sheet','ANOVA_PC1_PC2', 'Range','A2');

    writecell(TopLoadHeader, outXlsx, 'Sheet','TopLoadings', 'Range','A1');
    writecell(TopLoad,       outXlsx, 'Sheet','TopLoadings', 'Range','A2');

    fprintf('\n[OK] Done.\n');
    fprintf('[OK] Output folder: %s\n', outDir);
    fprintf('[OK] PCA summary: %s\n', outXlsx);
    fprintf('[OK] Log: %s\n', logFile);
    diary off;
end

% =========================================================================
% Helper functions (Objective 2)
% =========================================================================

function [DateCat, IrrCat, MngCat] = enforceCategoricals_obj2(T)
    % SamplingDate -> ordered categorical (chronological)
    d = T.SamplingDate;
    if ~isdatetime(d)
        d = datetime(d);
    end
    du = sort(unique(d));
    DateCat = categorical(d, du, cellstr(datestr(du,'yyyy-mm-dd')));

    % IrrigationRegime -> normalised labels (Irrigated, Rainfed) when possible
    raw = lower(strtrim(string(T.IrrigationRegime)));

    irrLabel = strings(size(raw));
    irrLabel(contains(raw,'irrig')) = "Irrigated";
    irrLabel(contains(raw,'rain'))  = "Rainfed";
    irrLabel(contains(raw,'riego')) = "Irrigated";
    irrLabel(contains(raw,'secano'))= "Rainfed";

    % Fallback: capitalise first letter for any unmatched label
    unmatched = (irrLabel=="");
    if any(unmatched)
        tmp = raw(unmatched);
        tmp = regexprep(tmp,'(^[a-z])','${upper($1)}');
        irrLabel(unmatched) = tmp;
    end

    IrrCat = categorical(irrLabel);
    % Force order if both present
    possibleCats = {'Irrigated','Rainfed'};
    existingCats = categories(IrrCat);
    toReorder = possibleCats(ismember(possibleCats, existingCats));
    if numel(toReorder)==2
        IrrCat = reordercats(IrrCat, toReorder);
    end

    % Vineyard floor management -> categorical (labels preserved)
    MngCat = categorical(string(T.VineyardFloorManagement));
end

function [X, wl] = extractSpectraPrefix(T, prefix)
    allVars = string(T.Properties.VariableNames);
    idx = startsWith(allVars, prefix);
    colNames = allVars(idx);
    if isempty(colNames)
        X = []; wl = [];
        return;
    end

    % Parse wavelength from variable names: FX10_nm_823p88 -> 823.88
    wlStr = replace(colNames, {prefix, 'p'}, {'','.'});
    wlVal = str2double(wlStr);
    [wl,ord] = sort(wlVal);
    colNames = colNames(ord);

    X = table2array(T(:, colNames));
end

function Xsnv = snvRows(X)
    mu = mean(X,2,'omitnan');
    sd = std(X,0,2,'omitnan');
    sd(sd==0) = 1;
    Xsnv = (X - mu) ./ sd;
end

function [coeff, score, explained] = pcaSVD(X, ncomp)
    % PCA via SVD on mean-centred matrix (robust for p >> n)
    mu = mean(X,1,'omitnan');
    Xc = X - mu;
    [U,S,V] = svd(Xc,'econ');
    r = min([ncomp, size(U,2), size(V,2)]);
    score = U(:,1:r) * S(1:r,1:r);
    coeff = V(:,1:r);

    latent = (diag(S).^2) ./ (size(Xc,1)-1);
    explainedAll = 100 * latent ./ sum(latent);
    explained = explainedAll(1:r);
end

function hf = plotScores_obj2(score, explained, groupVar, groupLabel, titleText)
    hf = figure('Color','w','Visible','off','Position',[100 100 900 650]);

    if isnumeric(groupVar); groupVar = string(groupVar); end
    if ~iscategorical(groupVar); groupVar = categorical(groupVar); end

    cats = categories(groupVar);
    C = lines(numel(cats));

    hold on;
    for i = 1:numel(cats)
        idx = groupVar == cats{i};
        scatter(score(idx,1), score(idx,2), 70, C(i,:), 'filled', ...
            'MarkerEdgeColor','k', 'MarkerFaceAlpha',0.85, 'DisplayName', string(cats{i}));
    end

    xlabel(sprintf('PC1 (%.1f%%)', explained(1)), 'FontName','Times New Roman','FontSize',12);
    ylabel(sprintf('PC2 (%.1f%%)', explained(2)), 'FontName','Times New Roman','FontSize',12);
    title({titleText; sprintf('PCA coloured by %s (mean of 70 berries per sample)', groupLabel)}, ...
        'FontName','Times New Roman','FontSize',14,'FontWeight','normal');

    lg = legend('Location','best');
    grid on; box on;
    set(gca,'FontName','Times New Roman','FontSize',10,'LineWidth',1.2,'Color','w','XColor','k','YColor','k');
    set(lg,'Box','on','EdgeColor','k','TextColor','k','Color','w','FontName','Times New Roman','FontSize',10);
    hold off;
end

function hf = plotLoadings_obj2(wl, coeff, explained, titleText)
    hf = figure('Color','w','Visible','off','Position',[100 100 900 450]);
    plot(wl, coeff(:,1), 'LineWidth', 1.8, 'DisplayName', sprintf('PC1 (%.1f%%)', explained(1)));
    hold on;
    plot(wl, coeff(:,2), 'LineWidth', 1.8, 'DisplayName', sprintf('PC2 (%.1f%%)', explained(2)));
    yline(0,'--k','LineWidth',1.0);
    xlabel('Wavelength (nm)','FontName','Times New Roman','FontSize',12);
    ylabel('Loading (a.u.)','FontName','Times New Roman','FontSize',12);
    title([titleText ' — loadings'], 'FontName','Times New Roman','FontSize',14,'FontWeight','normal');
    lg = legend('Location','best');
    grid on; box on;
    set(gca,'FontName','Times New Roman','FontSize',10,'LineWidth',1.2,'Color','w','XColor','k','YColor','k');
    set(lg,'Box','on','EdgeColor','k','TextColor','k','Color','w','FontName','Times New Roman','FontSize',10);
    hold off;
end

function hf = plotMeanAndDiffSpectra_obj2(wl, X, groupVar, titleText)
    hf = figure('Color','w','Visible','off','Position',[100 100 950 500]);
    if ~iscategorical(groupVar); groupVar = categorical(groupVar); end
    cats = categories(groupVar);
    C = lines(numel(cats));

    subplot(1,2,1);
    hold on;
    for i = 1:numel(cats)
        idx = groupVar == cats{i};
        mu = mean(X(idx,:),1,'omitnan');
        plot(wl, mu, 'LineWidth', 1.8, 'Color', C(i,:), 'DisplayName', cats{i});
    end
    xlabel('Wavelength (nm)','FontName','Times New Roman','FontSize',12);
    ylabel('Mean spectrum (a.u.)','FontName','Times New Roman','FontSize',12);
    title([titleText ' — group means'], 'FontName','Times New Roman','FontSize',13,'FontWeight','normal');
    grid on; box on;

    subplot(1,2,2);
    if numel(cats) == 2
        a = mean(X(groupVar==cats{1},:),1,'omitnan');
        b = mean(X(groupVar==cats{2},:),1,'omitnan');
        plot(wl, a-b, 'LineWidth', 1.8);
        yline(0,'--k','LineWidth',1.0);
        title(sprintf('%s — %s minus %s', titleText, cats{1}, cats{2}), ...
              'FontName','Times New Roman','FontSize',13,'FontWeight','normal');
        ylabel('Difference (a.u.)','FontName','Times New Roman','FontSize',12);
    else
        text(0.1,0.5,'Difference plot requires exactly 2 groups','FontName','Times New Roman','FontSize',12);
        axis off;
    end
    xlabel('Wavelength (nm)','FontName','Times New Roman','FontSize',12);
    grid on; box on;

    ax = findobj(hf,'Type','axes');
    set(ax,'FontName','Times New Roman','FontSize',10,'LineWidth',1.2,'Color','w','XColor','k','YColor','k');
end

function exportPub_obj2(fig, folder, name)
    % Export only .fig and .png (no .tif)
    savefig(fig, fullfile(folder, [name '.fig']));
    exportgraphics(fig, fullfile(folder, [name '.png']), 'Resolution', 300, 'BackgroundColor','white');
    close(fig);
end

function [p, eta2, ng] = anova1_p_eta(y, groupVar)
    if ~iscategorical(groupVar); groupVar = categorical(groupVar); end
    cats = categories(groupVar);
    ng = numel(cats);
    y = y(:);

    ok = isfinite(y) & ~isundefined(groupVar);
    y = y(ok);
    groupVar = groupVar(ok);
    cats = categories(groupVar);

    grand = mean(y);
    SST = sum((y - grand).^2);
    SSW = 0;
    nTot = numel(y);

    for i = 1:numel(cats)
        idx = groupVar == cats{i};
        yi = y(idx);
        if isempty(yi); continue; end
        mi = mean(yi);
        SSW = SSW + sum((yi - mi).^2);
    end

    SSB = SST - SSW;
    dfb = numel(cats) - 1;
    dfw = nTot - numel(cats);

    MSB = SSB / max(dfb,1);
    MSW = SSW / max(dfw,1);

    F = MSB / MSW;
    p = 1 - fcdf(F, dfb, dfw);
    eta2 = SSB / SST;
end
