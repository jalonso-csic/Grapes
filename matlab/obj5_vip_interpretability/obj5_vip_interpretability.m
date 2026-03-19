function obj5_vip_interpretability()
% obj5_vip_interpretability
% ------------------------------------------------------------
% Objective 5 / Results 3.5
% VIP-based interpretability for selected deployable traits.
%
% This workflow:
%   1. reads the master Excel matrix and data dictionary,
%   2. reconstructs the selected trait-specific model settings,
%   3. rebuilds predictor blocks for FX10, FX17, or FUSION,
%   4. applies the same preprocessing logic used in the supervised pipeline,
%   5. fits the final PLS-R model with the fixed selected LV,
%   6. computes VIP profiles,
%   7. extracts Top-5 VIP windows above the threshold,
%   8. exports VIP tables, per-trait figures, a multi-panel Figure 6,
%      and a run log.
%
% Output location:
%   <pwd>/Objetivo_5/
%     - tables/O5_VIP_Profiles_SelectedTraits.xlsx
%     - tables/Table6_VIP_Windows_Top5.xlsx
%     - figures/Figure6_VIP_Profiles_SelectedTraits.fig/.png
%     - figures/VIP_<EndpointVar>.fig/.png
%     - O5_VIP_runlog.txt
%
% Notes:
%   - Selected traits and model settings are fixed inside the script.
%   - White-background export is enforced for publication-safe figures.
%   - Single-point VIP windows are skipped during area-based ranking.
% ------------------------------------------------------------

rng(20260226,'twister');

% =================== Configuration ===================
cfg = struct();
cfg.inputFile  = fullfile(pwd,'Matriz_CHEM_HSI_MASTER_96.xlsx');
cfg.sheetData  = 'Matriz';
cfg.sheetDict  = 'DataDictionary';

cfg.outDir     = fullfile(pwd,'Objetivo_5');
cfg.figDir     = fullfile(cfg.outDir,'figures');
cfg.tabDir     = fullfile(cfg.outDir,'tables');
cfg.figureDPI  = 400;

% Preprocessing (match Objective 3)
cfg.sgPolyOrder = 2;
cfg.sgFrameLen  = 15;

% Fusion (match Objective 3)
cfg.fusionCut_nm = 950;
cfg.fusionL2norm = true;

% VIP windows
cfg.vipThreshold    = 1.0;
cfg.topKWindows     = 5;
cfg.primaryCumFrac  = 0.75;

% Shading for VIP windows (visual aid)
cfg.shadeTopWindows = true;
cfg.alphaPrimary    = 0.35;
cfg.alphaSecondary  = 0.18;
cfg.shadePrimaryColor   = [0.75 0.75 0.75];
cfg.shadeSecondaryColor = [0.90 0.90 0.90];

% ------------------- Deployable traits -------------------
vipTargets = table();
vipTargets.EndpointVar = { ...
    'TA_gL'
    'TSS_Brix'
    'TotalPhenolicPotential_mgkg'
    'MalicAcid_gL'
    'pH'
    'TotalAnthocyanins_mgkg'
    };

vipTargets.Sensor = { ...
    'FX10'
    'FUSION'
    'FX10'
    'FUSION'
    'FUSION'
    'FUSION'
    };

vipTargets.Preproc = { ...
    'SG2'
    'RAW'
    'SG2'
    'SNV'
    'SG1'
    'SG2'
    };

vipTargets.LV = [ ...
    5
    2
    7
    5
    3
    2
    ];

% Force as COLUMN variables (prevents "only 1 trait" bug)
vipTargets.EndpointVar = vipTargets.EndpointVar(:);
vipTargets.Sensor      = vipTargets.Sensor(:);
vipTargets.Preproc     = vipTargets.Preproc(:);
vipTargets.LV          = vipTargets.LV(:);

assert(height(vipTargets) == numel(vipTargets.LV), 'vipTargets size mismatch.');

% =================== I/O ===================
assert(isfile(cfg.inputFile), 'Input file not found: %s', cfg.inputFile);
T  = readtable(cfg.inputFile, 'Sheet', cfg.sheetData, 'PreserveVariableNames', true);
DD = readtable(cfg.inputFile, 'Sheet', cfg.sheetDict, 'PreserveVariableNames', true);

if ~exist(cfg.outDir,'dir'); mkdir(cfg.outDir); end
if ~exist(cfg.figDir,'dir'); mkdir(cfg.figDir); end
if ~exist(cfg.tabDir,'dir'); mkdir(cfg.tabDir); end

% Spectral columns
fx10Cols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'FX10_nm_'));
fx17Cols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames,'FX17_nm_'));
assert(~isempty(fx10Cols) && ~isempty(fx17Cols), 'FX10/FX17 spectral columns not found.');

fx10_nm = parseNmFromVarNames(fx10Cols, 'FX10_nm_');
fx17_nm = parseNmFromVarNames(fx17Cols, 'FX17_nm_');

% Add labels/units (if available)
vipTargets = addUnitsAndLabels(vipTargets, DD);

outProfilesXlsx = fullfile(cfg.tabDir, 'O5_VIP_Profiles_SelectedTraits.xlsx');
outWindowsXlsx  = fullfile(cfg.tabDir, 'Table6_VIP_Windows_Top5.xlsx');

% Clean outputs to avoid stale runs
if exist(outProfilesXlsx,'file'); delete(outProfilesXlsx); end
if exist(outWindowsXlsx,'file');  delete(outWindowsXlsx);  end

% =================== Main loop ===================
allWindows = table();
plotData   = struct([]);
kPlot = 0;

for i = 1:height(vipTargets)

    yVar   = vipTargets.EndpointVar{i};
    sensor = vipTargets.Sensor{i};
    prep   = vipTargets.Preproc{i};
    LV     = vipTargets.LV(i);
    traitLabel = vipTargets.Endpoint{i};

    assert(ismember(yVar, T.Properties.VariableNames), 'Endpoint not found: %s', yVar);

    yAll = T.(yVar);
    ok   = isfinite(yAll);
    if nnz(ok) < 10
        warning('Endpoint %s: too few non-missing values (%d). Skipping.', yVar, nnz(ok));
        continue;
    end

    [X, wl, blockLabel] = buildXandWavelengths(T, ok, sensor, fx10Cols, fx17Cols, fx10_nm, fx17_nm, cfg);
    y = yAll(ok);

    Xp = applyPreproc(X, wl, blockLabel, sensor, prep, cfg);
    Xz = zscoreLocal(Xp);

    LVfit = min([LV, size(Xz,1)-1, rank(Xz)]);
    if LVfit < 1
        warning('Endpoint %s: LVfit < 1. Skipping.', yVar);
        continue;
    end

    [~,~,XS,~,~,~,~,stats] = plsregress(Xz, y, LVfit); %#ok<ASGLU>
    vip = computeVIP_singleY(stats.W, XS, y);

    % Write VIP profile
    profTbl = table(wl(:), vip(:), blockLabel(:), 'VariableNames', {'wavelength_nm','VIP','Block'});
    writetable(profTbl, outProfilesXlsx, 'Sheet', safeSheetName(yVar));

    % VIP windows Top-5 (VIP≥1)
    winTbl = vipWindowsTopK(wl, vip, blockLabel, cfg.vipThreshold, cfg.topKWindows, cfg.primaryCumFrac);
    if ~isempty(winTbl)
        winTbl.Trait       = repmat(string(traitLabel), height(winTbl), 1);
        winTbl.EndpointVar = repmat(string(yVar), height(winTbl), 1);
        winTbl.BestModel   = repmat(string(sprintf('%s + %s (LV=%d)', sensor, prep, LVfit)), height(winTbl), 1);
        allWindows = [allWindows; winTbl]; %#ok<AGROW>
    end

    % Store for multi-panel plot
    kPlot = kPlot + 1;
    plotData(kPlot).TraitLabel = traitLabel; %#ok<AGROW>
    plotData(kPlot).EndpointVar = yVar; %#ok<AGROW>
    plotData(kPlot).wl         = wl; %#ok<AGROW>
    plotData(kPlot).vip        = vip; %#ok<AGROW>
    plotData(kPlot).block      = blockLabel; %#ok<AGROW>
    plotData(kPlot).winTbl     = winTbl; %#ok<AGROW>

    % Per-trait plot
    fig = figure('Color','w','Position',[100 100 1100 340]);
    plotVIPcurve(wl, vip, blockLabel, traitLabel, cfg.vipThreshold, winTbl, cfg);
    savefig(fig, fullfile(cfg.figDir, sprintf('VIP_%s.fig', yVar)));
    exportgraphics(fig, fullfile(cfg.figDir, sprintf('VIP_%s.png', yVar)), ...
        'Resolution', cfg.figureDPI, 'BackgroundColor','white');
    close(fig);

end

% Write windows (Table 6)
if ~isempty(allWindows)
    allWindows = sortrows(allWindows, {'Trait','Rank'});
    writetable(allWindows, outWindowsXlsx, 'Sheet', 'Table6');
end

% Figure 6: multi-panel VIP curves (2 x 3) with panel letters
fig = figure('Color','w','Position',[60 60 1250 720]);
tiledlayout(2,3,'Padding','compact','TileSpacing','compact');

panelLetters = {'A','B','C','D','E','F'};

for i = 1:min(6, numel(plotData))
    ax = nexttile; %#ok<NASGU>
    plotVIPcurve(plotData(i).wl, plotData(i).vip, plotData(i).block, plotData(i).TraitLabel, cfg.vipThreshold, plotData(i).winTbl, cfg);
    addPanelLetter(gca, panelLetters{i});
end

savefig(fig, fullfile(cfg.figDir, 'Figure6_VIP_Profiles_SelectedTraits.fig'));
exportgraphics(fig, fullfile(cfg.figDir, 'Figure6_VIP_Profiles_SelectedTraits.png'), ...
    'Resolution', cfg.figureDPI, 'BackgroundColor','white');
close(fig);

% Run log
logFile = fullfile(cfg.outDir, 'O5_VIP_runlog.txt');
fid = fopen(logFile,'w');
fprintf(fid, 'Objective 5 VIP interpretability run (v8)\n');
fprintf(fid, 'Timestamp: %s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'Input: %s\n', cfg.inputFile);
fprintf(fid, 'Traits requested: %d\n', height(vipTargets));
fprintf(fid, 'Traits processed: %d\n', numel(plotData));
fprintf(fid, 'VIP threshold: %.2f; TopK=%d; primaryCum=%.2f\n', cfg.vipThreshold, cfg.topKWindows, cfg.primaryCumFrac);
fprintf(fid, 'Outputs:\n  %s\n  %s\n  %s\n', outProfilesXlsx, outWindowsXlsx, fullfile(cfg.figDir,'Figure6_VIP_Profiles_SelectedTraits.png'));
fclose(fid);

disp('--- Objective 5 (VIP) complete ---');
disp(['Saved VIP profiles: ' outProfilesXlsx]);
disp(['Saved VIP windows:  ' outWindowsXlsx]);
disp(['Saved Figure 6:     ' fullfile(cfg.figDir,'Figure6_VIP_Profiles_SelectedTraits.png')]);

end

% =================== Helper functions ===================

function nm = parseNmFromVarNames(varNames, prefix)
nm = nan(numel(varNames),1);
for i=1:numel(varNames)
    s = erase(varNames{i}, prefix);
    s = strrep(s,'p','.');
    nm(i) = str2double(s);
end
end

function S = safeSheetName(s)
S = regexprep(s,'[^A-Za-z0-9_]','_');
if strlength(S) > 31
    S = extractBefore(S, 32);
end
end

function Xz = zscoreLocal(X)
mu = mean(X,1,'omitnan');
sd = std(X,0,1,'omitnan');
sd(sd==0) = 1;
Xz = (X - mu) ./ sd;
Xz(~isfinite(Xz)) = 0;
end

function vipTargets = addUnitsAndLabels(vipTargets, DD)
vipTargets.Endpoint = vipTargets.EndpointVar;
if ~all(ismember({'Variable','Unit (recommended)','Description (EN)'}, DD.Properties.VariableNames))
    return;
end
label = strings(height(vipTargets),1);
for i=1:height(vipTargets)
    v = vipTargets.EndpointVar{i};
    hit = strcmp(DD.Variable, v);
    if any(hit)
        u = string(DD.("Unit (recommended)")(find(hit,1,'first')));
        d = string(DD.("Description (EN)")(find(hit,1,'first')));
        if ismissing(u); u=""; end
        if ismissing(d); d=string(v); end
        if strlength(u)>0
            label(i) = d + " (" + u + ")";
        else
            label(i) = d;
        end
    else
        label(i) = string(v);
    end
end
vipTargets.Endpoint = cellstr(label);
end

function [X, wl, block] = buildXandWavelengths(T, ok, sensor, fx10Cols, fx17Cols, fx10_nm, fx17_nm, cfg)
sensor = upper(sensor);
switch sensor
    case 'FX10'
        X  = T{ok, fx10Cols};
        wl = fx10_nm(:);
        block = repmat("FX10", numel(wl), 1);

    case 'FX17'
        X  = T{ok, fx17Cols};
        wl = fx17_nm(:);
        block = repmat("FX17", numel(wl), 1);

    case 'FUSION'
        keep10 = fx10_nm <= cfg.fusionCut_nm;
        keep17 = fx17_nm >  cfg.fusionCut_nm;

        X10 = T{ok, fx10Cols(keep10)};
        X17 = T{ok, fx17Cols(keep17)};

        wl10 = fx10_nm(keep10);
        wl17 = fx17_nm(keep17);

        X  = [X10, X17];
        wl = [wl10(:); wl17(:)];
        block = [repmat("FX10", numel(wl10), 1); repmat("FX17", numel(wl17), 1)];
    otherwise
        error('Unknown sensor: %s', sensor);
end
end

function Xp = applyPreproc(X, wl, block, sensor, prep, cfg)
prep = upper(strrep(prep,' ', ''));
sensor = upper(sensor);

if strcmp(sensor,'FUSION')
    idx10 = (block=="FX10");
    idx17 = (block=="FX17");
    X10 = X(:, idx10);
    X17 = X(:, idx17);
    wl10 = wl(idx10);
    wl17 = wl(idx17);

    d10 = median(diff(wl10),'omitnan');
    d17 = median(diff(wl17),'omitnan');

    [X10p, X17p] = applyPreprocFusionBlocks(X10, X17, prep, cfg, d10, d17);

    if cfg.fusionL2norm
        X10p = l2normRows(X10p);
        X17p = l2normRows(X17p);
    end

    Xp = [X10p, X17p];
    return;
end

d = median(diff(wl),'omitnan');
switch prep
    case 'RAW'
        Xp = X;
    case 'SNV'
        Xp = snvRows(X);
    case 'SG1'
        Xp = sgolayDerivative2D(X, cfg.sgPolyOrder, cfg.sgFrameLen, 1, d);
    case 'SG2'
        Xp = sgolayDerivative2D(X, cfg.sgPolyOrder, cfg.sgFrameLen, 2, d);
    case 'SNV+SG1'
        Xp = sgolayDerivative2D(snvRows(X), cfg.sgPolyOrder, cfg.sgFrameLen, 1, d);
    case 'SNV+SG2'
        Xp = sgolayDerivative2D(snvRows(X), cfg.sgPolyOrder, cfg.sgFrameLen, 2, d);
    otherwise
        error('Unknown preprocessing: %s', prep);
end
end

function [X10p, X17p] = applyPreprocFusionBlocks(X10, X17, prep, cfg, d10, d17)
switch prep
    case 'RAW'
        X10p = X10; X17p = X17;
    case 'SNV'
        X10p = snvRows(X10); X17p = snvRows(X17);
    case 'SG1'
        X10p = sgolayDerivative2D(X10, cfg.sgPolyOrder, cfg.sgFrameLen, 1, d10);
        X17p = sgolayDerivative2D(X17, cfg.sgPolyOrder, cfg.sgFrameLen, 1, d17);
    case 'SG2'
        X10p = sgolayDerivative2D(X10, cfg.sgPolyOrder, cfg.sgFrameLen, 2, d10);
        X17p = sgolayDerivative2D(X17, cfg.sgPolyOrder, cfg.sgFrameLen, 2, d17);
    case 'SNV+SG1'
        X10p = sgolayDerivative2D(snvRows(X10), cfg.sgPolyOrder, cfg.sgFrameLen, 1, d10);
        X17p = sgolayDerivative2D(snvRows(X17), cfg.sgPolyOrder, cfg.sgFrameLen, 1, d17);
    case 'SNV+SG2'
        X10p = sgolayDerivative2D(snvRows(X10), cfg.sgPolyOrder, cfg.sgFrameLen, 2, d10);
        X17p = sgolayDerivative2D(snvRows(X17), cfg.sgPolyOrder, cfg.sgFrameLen, 2, d17);
    otherwise
        error('Unknown preprocessing: %s', prep);
end
end

function Xn = l2normRows(X)
den = sqrt(sum(X.^2,2));
den(den==0) = eps;
Xn = X ./ den;
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
filt = fliplr(g(:, derivOrder+1)');
left  = X(:, halfwin+1:-1:2);
right = X(:, end-1:-1:end-halfwin);
Xpad  = [left, X, right];
Xd = conv2(Xpad, filt, 'valid');
Xd = Xd * factorial(derivOrder) / (delta^derivOrder);
end

function vip = computeVIP_singleY(W, Tscore, y)
A = size(W,2);
p = size(W,1);
y = y(:);
y0 = y - mean(y,'omitnan');

SSY = zeros(A,1);
for a=1:A
    ta = Tscore(:,a);
    qa = (ta' * y0) / (ta' * ta);
    yhat_a = ta * qa;
    SSY(a) = sum(yhat_a.^2);
end

den = sum(SSY);
vip = zeros(p,1);
for j=1:p
    s = 0;
    for a=1:A
        wa = W(:,a);
        s = s + SSY(a) * (W(j,a)^2) / (wa' * wa);
    end
    vip(j) = sqrt(p * s / den);
end
end

function winTbl = vipWindowsTopK(wl, vip, block, thr, topK, primaryFrac)
blocks = unique(block,'stable');
rows = {};

for b = 1:numel(blocks)
    idx = (block==blocks(b));
    w = wl(idx);  w = w(:);
    v = vip(idx); v = v(:);

    mask = (v >= thr);
    if ~any(mask); continue; end

    d = diff([false; mask(:); false]);
    starts = find(d==1);
    ends   = find(d==-1)-1;

    for k=1:numel(starts)
        i1 = starts(k);
        i2 = ends(k);
        if i2 <= i1
            continue; % needs at least 2 points (avoid trapz overload)
        end
        w1 = w(i1); w2 = w(i2);
        area = trapz(w(i1:i2), v(i1:i2));
        vmean = mean(v(i1:i2));
        vmax  = max(v(i1:i2));
        rows = [rows; {string(blocks(b)), w1, w2, (w2-w1), area, vmean, vmax}]; %#ok<AGROW>
    end
end

if isempty(rows)
    winTbl = table();
    return;
end

winTbl = cell2table(rows, 'VariableNames', {'Block','nm_start','nm_end','width_nm','VIP_area','VIP_mean','VIP_max'});
winTbl = sortrows(winTbl, 'VIP_area', 'descend');

if height(winTbl) > topK
    winTbl = winTbl(1:topK,:);
end

winTbl.Rank = (1:height(winTbl))';
cum = cumsum(winTbl.VIP_area);
tot = sum(winTbl.VIP_area);
winTbl.Primary = cum <= (primaryFrac * tot);

winTbl = winTbl(:, {'Rank','Block','nm_start','nm_end','width_nm','VIP_area','VIP_mean','VIP_max','Primary'});
end

function plotVIPcurve(wl, vip, block, traitLabel, thr, winTbl, cfg, panelLetter)
if nargin < 8
    panelLetter = '';
end
if nargin < 7 || isempty(cfg)
    cfg = struct();
end
if ~isfield(cfg,'shadePrimaryOnly');  cfg.shadePrimaryOnly  = true; end
if ~isfield(cfg,'shadePrimaryAlpha'); cfg.shadePrimaryAlpha   = 0.45; end
if ~isfield(cfg,'shadeSecondaryAlpha'); cfg.shadeSecondaryAlpha = 0.00; end
if ~isfield(cfg,'shadePrimaryColor'); cfg.shadePrimaryColor   = [0.75 0.75 0.75]; end
if ~isfield(cfg,'shadeSecondaryColor'); cfg.shadeSecondaryColor = [0.90 0.90 0.90]; end
if nargin < 8
    panelLetter = '';
end
ax = gca;
if nargin < 6
    primaryOnly = false;
end
set(gcf,'Color','w');
set(ax,'Color','w','XColor','k','YColor','k','LineWidth',0.8);
hold on;

wl = wl(:);
vip = vip(:);

% Limits first (so shading uses fixed y-range)
xlim([min(wl) max(wl)]);
ymax = max([thr*1.2, max(vip)*1.05]);
ymax = ceil(ymax*10)/10;
ylim([0 ymax]);

% Shaded VIP windows (Top-5 VIP≥1)
if cfg.shadeTopWindows && ~isempty(winTbl)
    shadeVIPwindows(winTbl, cfg.shadePrimaryColor, cfg.shadeSecondaryColor, cfg.alphaPrimary, cfg.alphaSecondary);
end

% Plot VIP curves (by block)
blocks = unique(block,'stable');
for b = 1:numel(blocks)
    idx = (block==blocks(b));
    plot(wl(idx), vip(idx), 'LineWidth', 1.2);
end

% VIP=1 reference
yline(thr, '--', 'LineWidth', 1.0);

grid on;
ax.GridColor = [0.85 0.85 0.85];
ax.GridAlpha = 1;

xlabel('Wavelength (nm)');
ylabel('VIP');
% Panel letter (uppercase, no parentheses)
if ~isempty(panelLetter)
    text(0.02, 0.96, panelLetter, 'Units','normalized', 'FontWeight','bold', ...
        'FontName','Times New Roman', 'FontSize', 12, 'Color','k');
end


% Force integer x-axis start/end (no decimals)
xStart = 400;
if max(wl) > 1200
    xEnd = 1800;
else
    xEnd = 1000;
end
xlim([xStart xEnd]);
xticks(xStart:200:xEnd);

% Ensure y-limits defined before shading
ymax = max([thr*1.2, max(vip)*1.05]);
ymax = ceil(ymax*10)/10;
ylim([0 ymax]);

% Shade Primary VIP windows (Top-5, VIP>=threshold)
if ~isempty(winTbl) && istable(winTbl)
    shadeVIPwindows(winTbl, cfg.shadePrimaryColor, cfg.shadeSecondaryColor, cfg.shadePrimaryAlpha, cfg.shadeSecondaryAlpha, cfg.shadePrimaryOnly);
end




set(ax,'FontName','Times New Roman','FontSize',10);
set(get(ax,'XLabel'),'FontName','Times New Roman','FontSize',12);
set(get(ax,'YLabel'),'FontName','Times New Roman','FontSize',12);

end

function addPanelLetter(ax, letter)
% Adds (a)-(f) label at the top-left of the axes (normalized coords)
axes(ax); %#ok<LAXES>
text(0.02, 0.95, letter, 'Units','normalized', ...
    'FontName','Times New Roman', 'FontSize',12, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top', 'Color','k');
end


function shadeVIPwindows(winTbl, colPrimary, colSecondary, aPrimary, aSecondary, primaryOnly)
% Draw shaded rectangles for VIP windows.
% Expects winTbl columns: nm_start, nm_end, Primary (logical)
ax = gca;
if nargin < 6
    primaryOnly = false;
end
yl = ylim(ax);
for r = 1:height(winTbl)
    x1 = winTbl.nm_start(r);
    x2 = winTbl.nm_end(r);
    if ~isfinite(x1) || ~isfinite(x2) || x2<=x1
        continue;
    end
    if ismember('Primary', winTbl.Properties.VariableNames)
    isPrim = logical(winTbl.Primary(r));
else
    isPrim = true; % if no flag, treat as primary
end
if primaryOnly && ~isPrim
    continue;
end
if isPrim
    fc = colPrimary; fa = aPrimary;
else
    fc = colSecondary; fa = aSecondary;
end
    patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], fc, ...
        'FaceAlpha', fa, 'EdgeColor','none', 'HandleVisibility','off');
end
% Ensure patches stay behind lines
children = ax.Children;
% move patches to bottom
patchIdx = arrayfun(@(h) isa(h,'matlab.graphics.primitive.Patch'), children);
ax.Children = [children(~patchIdx); children(patchIdx)];
end