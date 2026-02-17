% BRC07_OBJ4_VIPStability_SelectedTiers_FigsAndBands_SNV_SG2nd_v1.m
% -------------------------------------------------------------------------
% Objective 4 — VIP stability mapping and operational band export
%
% This script consumes Objective 3 outputs (VIP raw stores + summary tables)
% and produces:
%   (1) Figure 1: VIP stability heatmap for Tier 1 + selected Tier 2 + selected Tier 3S
%       - Tier 1 shown at the top, then Tier 2, then Tier 3S
%       - Within each tier, endpoints are sorted by a chosen performance metric
%   (2) Figure 2: VIP stability line profiles for selected Tier 3S endpoints
%   (3) Excel export of stable wavelengths and stable contiguous regions
%
% Preprocessing context (fixed for this Objective 4 script):
%   SNV + Savitzky–Golay 2nd derivative (poly order = 2, frame length = 11)
%
% VIP stability definition:
%   VIP is normalised per model so that mean(VIP) across wavelengths = 1:
%       VIP_norm(model, :) = VIP(model, :) / mean(VIP(model, :))
%   Stability(λ) = proportion of models with VIP_norm(λ) > 1
%
% Stable region rule:
%   stability >= STAB_THR AND contiguous run length >= MIN_BANDS
%
% Requirements:
%   - MATLAB Statistics & Machine Learning Toolbox (for earlier Obj3 generation)
%   - Objective 3 outputs exist:
%       * OBJ3_VIPraw_Tier1_*SG2nd*.mat
%       * OBJ3_VIPraw_Tier2_*SG2nd*.mat
%       * OBJ3S_VIPraw_Tier3S_*SG2nd*.mat
%       * OBJ3_Summary_Tier1_*SG2nd*.xlsx (optional, for ordering)
%       * OBJ3_Summary_Tier2_*SG2nd*.xlsx (optional, for ordering)
%       * OBJ3S_Summary_Tier3S_*SG2nd*.xlsx (optional, for ordering)
%
% Notes:
%   - If summary workbooks are not found, the plotting order is kept as provided.
%   - Endpoint matching is robust to minor naming differences (e.g., '.' vs '_')
%     using matlab.lang.makeValidName().
% -------------------------------------------------------------------------

clear; clc;

%% =========================
% USER SETTINGS
% =========================

% Path to the Objective 3 output folder (the folder containing VIPraw .mat files).
% Leave "" to auto-detect, otherwise set an absolute/relative path.
OBJ3DIR_MANUAL = "";

% Output folder for Objective 4 products
OUTDIR = fullfile(pwd, 'OBJ4_OUT');
if ~exist(OUTDIR,'dir'); mkdir(OUTDIR); end

% Endpoints to include (edit as needed)
Tier1_All = ["Antihypertensive_act.", "ABTS", "DPPH", "Total_phenolics", "Extraction_yield"];

% Selected Tier 2 subset (edit as needed)
Tier2_Sel = ["SUM_GSL_indolic", "SUM_Amino_acids", "SUM_GSL_total"];

% Selected Tier 3S subset (edit as needed)
Tier3_Sel = ["Methoxyglucobrassicin_2", "Glucobrassicin", "Methoxyglucobrassicin_1"];

% Sorting metric for ordering endpoints within each tier (descending)
SORT_METRIC_COL = "R2_pooled";   % alternative: "RPD_pooled"

% Stability parameters
STAB_THR  = 0.70;   % stability threshold
MIN_BANDS = 25;     % minimum contiguous bands for a "region" (assumes 1 nm grid)

% Figure styling
FONT_MAIN = 'Times New Roman';
FS_AX     = 10;
FS_LAB    = 12;
FS_TTL    = 12;
PNG_DPI   = 400;

%% =========================
% LOCATE OBJ3 OUTPUT FOLDER
% =========================
OBJ3DIR = string(OBJ3DIR_MANUAL);

if strlength(OBJ3DIR)==0 || ~(exist(OBJ3DIR,'dir')==7)

    % Try common relative locations from current folder
    candDirs = { ...
        fullfile(pwd, 'OBJ3_OUT'); ...
        fullfile(pwd, '..', 'OBJ3_OUT'); ...
        fullfile(pwd, 'OBJ3', 'OBJ3_OUT'); ...
        fullfile(pwd, '..', 'OBJ3', 'OBJ3_OUT'); ...
        fullfile(pwd, 'OBJ3_OUT_4CFG'); ...
        fullfile(pwd, '..', 'OBJ3_OUT_4CFG') ...
    };

    found = false;
    for i=1:numel(candDirs)
        if exist(candDirs{i},'dir')==7
            OBJ3DIR = string(candDirs{i});
            found = true;
            break;
        end
    end

    if ~found
        fprintf('\nObjective 3 output folder was not found automatically.\n');
        tmp = uigetdir(pwd, 'Select your Objective 3 output folder (contains VIPraw .mat files)');
        assert(~isequal(tmp,0), 'No Objective 3 output folder selected.');
        OBJ3DIR = string(tmp);
    end
end

fprintf('OBJ3DIR = %s\n', OBJ3DIR);

%% =========================
% ORDER ENDPOINTS (OPTIONAL; REQUIRES SUMMARY XLSX)
% =========================
SUM_T1 = pickSummary(char(OBJ3DIR), '*OBJ3_Summary_Tier1*SG2nd*.xlsx',   'Tier 1 summary (SNV+SG2nd)');
SUM_T2 = pickSummary(char(OBJ3DIR), '*OBJ3_Summary_Tier2*SG2nd*.xlsx',   'Tier 2 summary (SNV+SG2nd)');
SUM_T3 = pickSummary(char(OBJ3DIR), '*OBJ3S_Summary_Tier3S*SG2nd*.xlsx', 'Tier 3S summary (SNV+SG2nd)');

Tier1_All = orderByMetric(Tier1_All, SUM_T1, SORT_METRIC_COL);
Tier2_Sel = orderByMetric(Tier2_Sel, SUM_T2, SORT_METRIC_COL);
Tier3_Sel = orderByMetric(Tier3_Sel, SUM_T3, SORT_METRIC_COL);

% Final plotting order: Tier 1 (top) -> Tier 2 -> Tier 3S (bottom)
labelsRaw = [Tier1_All, Tier2_Sel, Tier3_Sel];
nY = numel(labelsRaw);

%% =========================
% LOAD VIPraw STORES (SNV+SG2nd)
% =========================
MAT_T1 = pickMat(char(OBJ3DIR), '*Tier1*SNV*SG2nd*.mat',   'Tier 1 VIPraw (SNV+SG2nd)');
MAT_T2 = pickMat(char(OBJ3DIR), '*Tier2*SNV*SG2nd*.mat',   'Tier 2 VIPraw (SNV+SG2nd)');
MAT_T3 = pickMat(char(OBJ3DIR), '*Tier3S*SNV*SG2nd*.mat',  'Tier 3S VIPraw (SNV+SG2nd)');

VIP1 = loadVipStore(MAT_T1);
VIP2 = loadVipStore(MAT_T2);
VIP3 = loadVipStore(MAT_T3);

%% =========================
% BUILD STABILITY MATRIX
% =========================
[wl, p] = getWavelengthFromAny(VIP1, VIP2, VIP3, labelsRaw);
StabMat = nan(nY, p);

for i = 1:nY
    y = labelsRaw(i);
    [vip, wl_local] = fetchVIP(VIP1, VIP2, VIP3, y);
    assert(numel(wl_local)==p, 'Wavelength length mismatch for endpoint: %s', y);

    stab = vipStability(vip);   % returns [p x 1]
    StabMat(i,:) = stab(:)';
end

%% =========================
% FIGURE 1 — STABILITY HEATMAP
% =========================
fig1 = figure('Color','w','Units','pixels','Position',[80 80 1450 640]);
imagesc(wl, 1:nY, StabMat);
set(gca,'YDir','reverse');  % row 1 at TOP
caxis([0 1]);

cb = colorbar;
ylabel(cb, 'VIP stability', 'FontName',FONT_MAIN,'FontSize',FS_AX);

yticks(1:nY);
yticklabels(prettyLabels(labelsRaw));
set(gca,'FontName',FONT_MAIN,'FontSize',FS_AX,'TickLabelInterpreter','none');

xlabel('Wavelength (nm)', 'FontName',FONT_MAIN,'FontSize',FS_LAB);
title('VIP stability map (SNV + Savitzky–Golay 2nd derivative) — Tier 1 + selected Tier 2 + selected Tier 3S', ...
    'FontName',FONT_MAIN,'FontSize',FS_TTL);

drawnow;
saveas(fig1, fullfile(OUTDIR, 'Fig_OBJ4_VIPStability_Heatmap_T1_T2sel_T3sel_SNV_SG2nd.fig'));
exportgraphics(fig1, fullfile(OUTDIR, 'Fig_OBJ4_VIPStability_Heatmap_T1_T2sel_T3sel_SNV_SG2nd.png'), 'Resolution', PNG_DPI);

%% =========================
% FIGURE 2 — STABILITY PROFILES (TIER 3S SELECTED)
% =========================
fig2 = figure('Color','w','Units','pixels','Position',[80 80 1450 520]);
hold on;

for i = 1:numel(Tier3_Sel)
    y = Tier3_Sel(i);
    [vip, wl_local] = fetchVIP(VIP1, VIP2, VIP3, y);
    stab = vipStability(vip);
    plot(wl_local, stab, 'LineWidth', 1.4, 'DisplayName', prettyLabels(y));
end

yline(STAB_THR, '--', 'LineWidth', 1.0, ...
    'DisplayName', sprintf('Stability threshold (%.2f)', STAB_THR));

ylim([0 1.05]);
xlim([wl(1) wl(end)]);

set(gca,'FontName',FONT_MAIN,'FontSize',FS_AX);
xlabel('Wavelength (nm)', 'FontName',FONT_MAIN,'FontSize',FS_LAB);
ylabel('VIP stability', 'FontName',FONT_MAIN,'FontSize',FS_LAB);
title('VIP stability profiles (SNV + Savitzky–Golay 2nd derivative) — selected Tier 3S endpoints', ...
    'FontName',FONT_MAIN,'FontSize',FS_TTL);

lgd = legend('Location','southoutside','NumColumns',2);
set(lgd,'FontName',FONT_MAIN,'FontSize',FS_AX);
box on; grid off;

drawnow;
saveas(fig2, fullfile(OUTDIR, 'Fig_OBJ4_VIPStability_Lines_Tier3Sel_SNV_SG2nd.fig'));
exportgraphics(fig2, fullfile(OUTDIR, 'Fig_OBJ4_VIPStability_Lines_Tier3Sel_SNV_SG2nd.png'), 'Resolution', PNG_DPI);

%% =========================
% EXCEL EXPORT — STABLE BANDS + STABLE REGIONS + COUNTS
% =========================
stableBands   = table();
stableRegions = table();

counts = table('Size',[nY 4], ...
    'VariableTypes', {'string','double','double','double'}, ...
    'VariableNames', {'Endpoint','nStableNm','nRegions','maxStability'});

for i = 1:nY
    y = labelsRaw(i);
    yLab = string(prettyLabels(y));

    [vip, wl_local] = fetchVIP(VIP1, VIP2, VIP3, y);
    [stab, meanVIP] = vipStability_full(vip);  % includes meanVIP for reporting

    % (A) Stable wavelengths
    idx = find(stab >= STAB_THR);
    if ~isempty(idx)
        Tb = table();
        Tb.Endpoint  = repmat(yLab, numel(idx), 1);
        Tb.nm        = wl_local(idx);
        Tb.stability = stab(idx);
        Tb.meanVIP   = meanVIP(idx);
        stableBands  = [stableBands; Tb]; %#ok<AGROW>
    end

    % (B) Stable regions (contiguous runs)
    Tr = regionsFromMask(wl_local, stab, meanVIP, STAB_THR, MIN_BANDS);
    if ~isempty(Tr)
        Tr.Endpoint = repmat(yLab, height(Tr), 1);
        stableRegions = [stableRegions; Tr]; %#ok<AGROW>
    end

    counts.Endpoint(i)     = yLab;
    counts.nStableNm(i)    = numel(idx);
    counts.nRegions(i)     = height(Tr);
    counts.maxStability(i) = max(stab, [], 'omitnan');
end

xlsxOut = fullfile(OUTDIR, 'OBJ4_StableBands_SelectedEndpoints_SNV_SG2nd.xlsx');
writetable(stableRegions, xlsxOut, 'Sheet', 'StableRegions');
writetable(stableBands,   xlsxOut, 'Sheet', 'StableBands_nm');
writetable(counts,        xlsxOut, 'Sheet', 'Counts');

fprintf('\nObjective 4 outputs written to:\n  %s\n', OUTDIR);
fprintf('Excel stable bands written to:\n  %s\n\n', xlsxOut);

%% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function matPath = pickMat(folder, pattern, label)
    cand = dir(fullfile(folder, pattern));
    if isempty(cand)
        fprintf('\nCould not auto-detect %s MAT in:\n%s\n', label, folder);
        [f,p] = uigetfile('*.mat', ['Select ' label]);
        assert(~isequal(f,0), 'No MAT selected for %s.', label);
        matPath = fullfile(p,f);
    else
        matPath = fullfile(cand(1).folder, cand(1).name);
    end
    fprintf('Using %s: %s\n', label, matPath);
end

function xlsxPath = pickSummary(folder, pattern, label)
    cand = dir(fullfile(folder, pattern));
    if isempty(cand)
        warning('%s not found (pattern "%s" in %s). Keeping the original order for that tier.', ...
            label, pattern, folder);
        xlsxPath = "";
        return;
    end
    xlsxPath = fullfile(cand(1).folder, cand(1).name);
    fprintf('Using %s: %s\n', label, xlsxPath);
end

function yOut = orderByMetric(yIn, summaryPath, metricCol)
    yIn = string(yIn);
    if strlength(string(summaryPath))==0 || exist(summaryPath,'file')~=2
        yOut = yIn;
        return;
    end

    S = readtable(summaryPath);
    if ~ismember('Y', S.Properties.VariableNames) || ~ismember(metricCol, S.Properties.VariableNames)
        warning('Summary file missing expected columns (Y, %s): %s. Keeping original order.', metricCol, summaryPath);
        yOut = yIn;
        return;
    end

    SY  = string(S.Y);
    SYv = string(arrayfun(@(s) matlab.lang.makeValidName(s), SY, 'UniformOutput', false));

    score = nan(size(yIn));

    for i=1:numel(yIn)
        q = yIn(i);

        % 1) exact match
        idx = strcmp(SY, q);
        if any(idx)
            score(i) = S.(metricCol)(find(idx,1,'first'));
            continue;
        end

        % 2) match by makeValidName equivalence (handles '.' vs '_' etc.)
        qv = string(matlab.lang.makeValidName(q));
        idx2 = strcmp(SYv, qv);
        if any(idx2)
            score(i) = S.(metricCol)(find(idx2,1,'first'));
            continue;
        end
    end

    scoreTmp = score;
    scoreTmp(isnan(scoreTmp)) = -inf;      % not found -> last
    [~, ord] = sort(scoreTmp, 'descend');
    yOut = yIn(ord);
end

function VIP_STORE = loadVipStore(matPath)
    S = load(matPath);
    assert(isfield(S,'VIP_STORE'), 'MAT file does not contain VIP_STORE: %s', matPath);
    VIP_STORE = S.VIP_STORE;
end

function [wl, p] = getWavelengthFromAny(VIP1, VIP2, VIP3, labelsRaw)
    wl = [];
    for i=1:numel(labelsRaw)
        yv = matlab.lang.makeValidName(labelsRaw(i));
        if isfield(VIP1, yv); wl = double(VIP1.(yv).wl(:)); break; end
        if isfield(VIP2, yv); wl = double(VIP2.(yv).wl(:)); break; end
        if isfield(VIP3, yv); wl = double(VIP3.(yv).wl(:)); break; end
    end
    assert(~isempty(wl), 'Could not find wavelength vector (wl) in any VIP_STORE.');
    p = numel(wl);
end

function [vip, wl] = fetchVIP(VIP1, VIP2, VIP3, y)
    yv = matlab.lang.makeValidName(y);

    if isfield(VIP1, yv)
        vip = double(VIP1.(yv).vip);
        wl  = double(VIP1.(yv).wl(:));
        vip = orientVIP(vip, numel(wl));
        return;
    end
    if isfield(VIP2, yv)
        vip = double(VIP2.(yv).vip);
        wl  = double(VIP2.(yv).wl(:));
        vip = orientVIP(vip, numel(wl));
        return;
    end
    if isfield(VIP3, yv)
        vip = double(VIP3.(yv).vip);
        wl  = double(VIP3.(yv).wl(:));
        vip = orientVIP(vip, numel(wl));
        return;
    end

    error('Endpoint not found in VIP stores: %s', y);
end

function vip = orientVIP(vip, p)
    % Expected orientation: [nModels x p]
    % If stored as [p x nModels], transpose.
    if size(vip,2) ~= p && size(vip,1) == p
        vip = vip.';
    end
    assert(size(vip,2)==p, 'VIP has unexpected size; expected p=%d columns.', p);
end

function stab = vipStability(vip)
    % Returns stability vector [p x 1]
    mu = mean(vip, 2, 'omitnan');     % per-model mean across wavelengths
    mu(mu==0) = 1;
    vipN = vip ./ mu;                  % normalised VIP (mean=1 per model)
    stab = mean(vipN > 1, 1, 'omitnan');% proportion across models
    stab = stab(:);
end

function [stab, meanVIP] = vipStability_full(vip)
    % Returns stability [p x 1] and meanVIP [p x 1] for reporting
    mu = mean(vip, 2, 'omitnan');
    mu(mu==0) = 1;
    vipN = vip ./ mu;
    stab = mean(vipN > 1, 1, 'omitnan');
    meanVIP = mean(vipN, 1, 'omitnan');
    stab = stab(:);
    meanVIP = meanVIP(:);
end

function Tr = regionsFromMask(wl, stab, meanVIP, thr, minBands)
    mask = stab >= thr;
    Tr = table();
    if ~any(mask); return; end

    idx = find(mask);
    starts = idx([true; diff(idx) > 1]);
    ends   = idx([diff(idx) > 1; true]);

    keep = (ends - starts + 1) >= minBands;
    starts = starts(keep);
    ends   = ends(keep);

    if isempty(starts); return; end

    nR = numel(starts);
    Tr = table('Size',[nR 8], ...
        'VariableTypes', {'double','double','double','double','double','double','double','double'}, ...
        'VariableNames', {'start_nm','end_nm','width_nm','n_bands','mean_stability','max_stability','mean_VIP','max_meanVIP'});

    for r = 1:nR
        s = starts(r); e = ends(r);
        seg = s:e;

        Tr.start_nm(r) = wl(s);
        Tr.end_nm(r)   = wl(e);
        Tr.width_nm(r) = wl(e) - wl(s);
        Tr.n_bands(r)  = numel(seg);

        Tr.mean_stability(r) = mean(stab(seg), 'omitnan');
        Tr.max_stability(r)  = max(stab(seg), [], 'omitnan');
        Tr.mean_VIP(r)       = mean(meanVIP(seg), 'omitnan');
        Tr.max_meanVIP(r)    = max(meanVIP(seg), [], 'omitnan');
    end
end

function out = prettyLabels(names)
    names = string(names);
    out = strings(size(names));

    for i=1:numel(names)
        s = names(i);

        % Labels only: remove "SUM_"
        s = regexprep(s, '^(?i)SUM_', '');

        % Replace underscores with spaces
        s = replace(s, "_", " ");
        s = strtrim(regexprep(s, '\s+', ' '));

        % Targeted pretty names
        if contains(lower(s), "antihypertensive")
            s = "Antihypertensive activity";
        end
        if strcmpi(s, "GSL indolic"); s = "Indolic glucosinolates"; end
        if strcmpi(s, "GSL total");   s = "Total glucosinolates"; end

        out(i) = s;
    end
end
