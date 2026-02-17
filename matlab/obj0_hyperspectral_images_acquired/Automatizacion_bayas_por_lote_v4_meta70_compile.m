function Automatizacion_bayas_por_lote_v4_meta70_compile()
% Automatizacion_bayas_por_lote_v4_meta70_compile
%
% FX17 (SWIR) batch pipeline for berry-level spectral extraction.
%
% Pipeline (algorithmic intent preserved):
%  1) Reflectance calibration using sample-specific WHITE/DARK references.
%  2) HSV pre-mask ("anti-rachis") to suppress rachis/background.
%  3) SAM segmentation on pseudo-RGB; restrict labels to HSV mask.
%  4) Split oversized blobs via watershed; re-label berries.
%  5) Select top 90 berries by area; ROI = full berry eroded by 2 px (fallback to full berry).
%  6) Spectral QC selection: keep 70 closest to the median spectrum (Euclidean distance).
%  7) Export per-sample Excel tables + QC figures, then compile across samples.
%
% Expected structure (run from experiment root):
%   ROOT/
%     FX17_.../
%       capture/
%         <sample>.hdr/.raw
%         WHITEREF_<sample>.hdr/.raw (or WHITEREF_FX17_<sample>.*)
%         DARKREF_<sample>.hdr/.raw  (or DARKREF_FX17_<sample>.*)
%
% Outputs (created in ROOT):
%   - results_STEP1_2_FX17_AntiRachis_SAM/
%   - results_Spectra90_FX17_CentralROIs/
%   - results_Spectra70_FX17_CentralROIs/
%
% Requirements:
%   - MATLAB + Image Processing Toolbox
%   - Add-on: Image Processing Toolbox Model for Segment Anything Model (imsegsam)
%   - ENVI readers: read_envi_info/read_envi_cube OR enviinfo_local/enviread_local
%
close all; clc;

rootDir = pwd;

% --- SAM availability ---
if ~exist('imsegsam','file')
    error(['Missing imsegsam. Install the add-on "Image Processing Toolbox Model ' ...
           'for Segment Anything Model" from Add-Ons.']);
end

% --- Bind ENVI readers ---
read_info = @read_envi_info;
read_cube = @read_envi_cube;
if exist('enviinfo_local','file'), read_info = @enviinfo_local; end
if exist('enviread_local','file'), read_cube = @enviread_local; end

% --- Output folders ---
outStep12Dir    = fullfile(rootDir, 'results_STEP1_2_FX17_AntiRachis_SAM');
outSpectra90Dir = fullfile(rootDir, 'results_Spectra90_FX17_CentralROIs');
outSpectra70Dir = fullfile(rootDir, 'results_Spectra70_FX17_CentralROIs');
ensureDir(outStep12Dir);
ensureDir(outSpectra90Dir);
ensureDir(outSpectra70Dir);

% --- Discover samples ---
d = dir(rootDir);
isSample = @(x) x.isdir && startsWith(x.name,'FX17_') && ~contains(lower(x.name),'resultado');
sampleDirs = d(arrayfun(isSample,d));
if isempty(sampleDirs)
    error('No FX17_* sample folders found in: %s', rootDir);
end
fprintf('FX17 samples found: %d\n', numel(sampleDirs));

% ---------- Parameters ----------
MIN_BERRY_AREA = 50;

HUE_RANGE = [0.55 0.8];
MIN_SATURATION_THRESHOLD = 0.35;
MAX_VALUE_THRESHOLD      = 0.85;

ERODE_PIX = 2;

design = makeDesignTable();

% ========= Main loop =========
for k = 1:numel(sampleDirs)
    sampleName = sampleDirs(k).name;
    captureDir = fullfile(rootDir, sampleName, 'capture');

    fprintf('\n=====================================================\n');
    fprintf('Processing FX17 sample: %s (%d/%d)\n', sampleName, k, numel(sampleDirs));

    % --- Paths ---
    hdr_sample = fullfile(captureDir, [sampleName '.hdr']);
    raw_sample = fullfile(captureDir, [sampleName '.raw']);

    % WHITE (two naming conventions)
    hdr_white = fullfile(captureDir, ['WHITEREF_' sampleName '.hdr']);
    raw_white = fullfile(captureDir, ['WHITEREF_' sampleName '.raw']);
    if ~exist(hdr_white,'file')
        hdr_white = fullfile(captureDir, ['WHITEREF_FX17_' sampleName '.hdr']);
        raw_white = fullfile(captureDir, ['WHITEREF_FX17_' sampleName '.raw']);
    end

    % DARK (two naming conventions)
    hdr_dark = fullfile(captureDir, ['DARKREF_' sampleName '.hdr']);
    raw_dark = fullfile(captureDir, ['DARKREF_' sampleName '.raw']);
    if ~exist(hdr_dark,'file')
        hdr_dark = fullfile(captureDir, ['DARKREF_FX17_' sampleName '.hdr']);
        raw_dark = fullfile(captureDir, ['DARKREF_FX17_' sampleName '.raw']);
    end

    if ~exist(hdr_white,'file') || ~exist(hdr_dark,'file')
        fprintf('  ERROR: WHITE/DARK references not found for %s\n', sampleName);
        continue;
    end

    % --- Load + calibrate ---
    try
        info  = read_info(hdr_sample);
        cube  = read_cube(info, raw_sample);

        infoW = read_info(hdr_white);
        white = read_cube(infoW, raw_white);

        infoD = read_info(hdr_dark);
        dark  = read_cube(infoD, raw_dark);

        [cube_cal, nb] = calibrateReflectance(cube, white, dark);
    catch ME
        fprintf('  ERROR (load/calibrate): %s\n', ME.message);
        continue;
    end

    % --- Wavelength vector ---
    wl = getWavelengthsFromInfoOrHdr(info, hdr_sample, nb);

    % --- Pseudo-RGB from calibrated cube (fixed relative positions) ---
    idxR = max(1, round(0.85*nb));
    idxG = max(1, round(0.55*nb));
    idxB = max(1, round(0.25*nb));
    rgb_pseudo = cat(3, cube_cal(:,:,idxR), cube_cal(:,:,idxG), cube_cal(:,:,idxB));
    rgb_pseudo = rgb_pseudo - min(rgb_pseudo(:));
    mx = max(rgb_pseudo(:));
    if mx>0, rgb_pseudo = rgb_pseudo./mx; end
    rgb_pseudo = imadjust(rgb_pseudo, stretchlim(rgb_pseudo,[0.01 0.99]));

    % =======================================================
    % STEP 1: HSV pre-mask
    % =======================================================
    fprintf('  STEP 1: HSV pre-mask\n');

    hsv = rgb2hsv(rgb_pseudo);
    hCh = hsv(:,:,1); sCh = hsv(:,:,2); vCh = hsv(:,:,3);

    mask_hue = (hCh > HUE_RANGE(1)) & (hCh < HUE_RANGE(2));

    tS = graythresh(sCh);
    mask_sat = (sCh > tS) & (sCh > MIN_SATURATION_THRESHOLD);

    tV = graythresh(vCh);
    mask_val_low  = (vCh > tV*0.5);
    mask_val_high = (vCh < MAX_VALUE_THRESHOLD);

    berryMask = mask_hue & mask_sat & mask_val_low & mask_val_high;

    berryMask = bwareaopen(berryMask, MIN_BERRY_AREA);
    berryMask = imclose(berryMask, strel('disk',3));
    berryMask = imfill(berryMask,'holes');
    berryMask = imopen(berryMask, strel('disk',1));
    berryMask = imdilate(berryMask, strel('disk',1));

    h = figure('Visible','off','Color','w'); imshow(berryMask);
    title(sprintf('STEP 1: HSV pre-mask — %s', sampleName), 'Interpreter','none');
    exportgraphics(h, fullfile(outStep12Dir, sprintf('STEP1_HSV_PreMask_%s.png', sampleName)), 'Resolution', 300);
    close(h);

    % =======================================================
    % STEP 2: SAM segmentation restricted to HSV mask + watershed split
    % =======================================================
    fprintf('  STEP 2: SAM segmentation (restricted)\n');

    I_sam = im2uint8(rgb_pseudo);
    masksSAM = imsegsam(I_sam);
    Lsam_raw = labelmatrix(masksSAM);

    % Restrict labels to HSV mask
    Lsam = Lsam_raw;
    Lsam(~berryMask) = 0;

    labels_inside = unique(Lsam(:));
    labels_inside(labels_inside==0) = [];
    if isempty(labels_inside)
        fprintf('  WARNING: No SAM regions inside HSV mask for %s\n', sampleName);
        continue;
    end

    areas_raw = zeros(numel(labels_inside),1);
    for i=1:numel(labels_inside)
        areas_raw(i) = nnz(Lsam==labels_inside(i));
    end
    medianArea = median(areas_raw);
    MIN_AREA_BERRY  = 0.25 * medianArea;
    MAX_AREA_SIMPLE = 1.8  * medianArea;

    [rows, cols, ~] = size(I_sam);
    Lberries = zeros(rows, cols, 'uint16');
    new_id = 0;

    for i=1:numel(labels_inside)
        lab = labels_inside(i);
        mask_lab = (Lsam==lab);
        A = nnz(mask_lab);
        if A==0, continue; end

        if A <= MAX_AREA_SIMPLE
            new_id = new_id + 1;
            Lberries(mask_lab) = new_id;
        else
            D = bwdist(~mask_lab);
            D = -D; D(~mask_lab) = Inf;
            Lw = watershed(D);

            subLabs = unique(Lw(mask_lab));
            subLabs(subLabs==0) = [];

            for s=1:numel(subLabs)
                ms = mask_lab & (Lw==subLabs(s));
                As = nnz(ms);
                if As < MIN_AREA_BERRY, continue; end
                new_id = new_id + 1;
                Lberries(ms) = new_id;
            end
        end
    end

    numBerries = new_id;
    if numBerries==0
        fprintf('  WARNING: No berries after refinement for %s\n', sampleName);
        continue;
    end

    h = figure('Visible','off','Color','w');
    imshow(label2rgb(Lberries,'jet','w','shuffle'));
    title(sprintf('STEP 2: %d berries (refined) — %s', numBerries, sampleName), 'Interpreter','none');
    exportgraphics(h, fullfile(outStep12Dir, sprintf('STEP2_Berries_Refined_%s.png', sampleName)), 'Resolution', 300);
    close(h);

    save(fullfile(outStep12Dir, sprintf('labels_SAM_%s.mat', sampleName)), 'Lberries');

    % White background outside HSV mask for clean overlays
    rgb_clean = rgb_pseudo;
    for c=1:3
        ch = rgb_clean(:,:,c);
        ch(~berryMask) = 1;
        rgb_clean(:,:,c) = ch;
    end

    % =======================================================
    % STEP 3: Select top 90 by area + ROI = eroded berry (2 px)
    % =======================================================
    fprintf('  STEP 3: select top 90 by area + extract spectra\n');

    stats = regionprops(Lberries,'Area');
    Area_px = reshape([stats.Area], [], 1);
    [~, idxSort] = sort(Area_px,'descend');
    Nkeep = min(90, numel(idxSort));
    idxTop = idxSort(1:Nkeep);

    seErode = strel('disk', ERODE_PIX);

    L_roi = zeros(rows, cols, 'uint16'); % ROI labels 1..Nkeep
    Spectra90 = zeros(Nkeep, nb);

    cube2 = reshape(cube_cal, [], nb);
    for ii=1:Nkeep
        idBerry = idxTop(ii);
        berry_mask = (Lberries==idBerry);

        roi_mask = imerode(berry_mask, seErode);
        if nnz(roi_mask) < 20
            roi_mask = berry_mask; % fallback
        end

        L_roi(roi_mask) = ii;
        roi_idx = roi_mask(:);
        Spectra90(ii,:) = mean(cube2(roi_idx,:),1);
    end

    % =======================================================
    % STEP 4: Export 90 + select/export 70 + QC figures
    % =======================================================
    fprintf('  STEP 4: export 90 + select/export 70\n');

    medSpec = median(Spectra90,1);
    dist_med = sqrt(sum((Spectra90 - medSpec).^2,2));
    [~, idxByDist] = sort(dist_med,'ascend');
    Nclean = min(70, Nkeep);
    idxClean = idxByDist(1:Nclean);

    wlLabels = makeNmLabels_dec2(wl);

    % 90
    T90 = table(repmat(string(sampleName), Nkeep,1), (1:Nkeep).', dist_med(:), ...
        'VariableNames', {'SampleName','ROI_ID','Dist_to_Median'});
    for j=1:nb
        T90.(wlLabels{j}) = Spectra90(:,j);
    end
    writetable(T90, fullfile(outSpectra90Dir, sprintf('Spectra90_CentralROIs_%s.xlsx', sampleName)), ...
        'Sheet','Spectra90','WriteMode','overwritesheet');

    % 70 with metadata
    Spectra70 = Spectra90(idxClean,:);
    dist70 = dist_med(idxClean);

    [Sensor, Variety, DateStr, SampleIndex, Treatment, Block, IrrigationRegime, VineyardFloorManagement] = ...
        parseMetaFromSampleName(sampleName, design);

    n70 = numel(idxClean);
    T70 = table( ...
        repmat(string(Sensor), n70,1), ...
        repmat(string(Variety), n70,1), ...
        repmat(string(DateStr), n70,1), ...
        repmat(double(SampleIndex), n70,1), ...
        repmat(string(Treatment), n70,1), ...
        repmat(string(Block), n70,1), ...
        repmat(string(IrrigationRegime), n70,1), ...
        repmat(string(VineyardFloorManagement), n70,1), ...
        repmat(string(sampleName), n70,1), ...
        idxClean(:), dist70(:), ...
        'VariableNames', {'Sensor','Variety','SamplingDate','SampleIndex','Treatment','Block','IrrigationRegime','VineyardFloorManagement', ...
                          'SampleName','ROI_ID','Dist_to_Median'} );
    for j=1:nb
        T70.(wlLabels{j}) = Spectra70(:,j);
    end
    writetable(T70, fullfile(outSpectra70Dir, sprintf('Spectra70_CentralROIs_%s.xlsx', sampleName)), ...
        'Sheet','Spectra70','WriteMode','overwritesheet');

    % QC figures
    h = figure('Visible','off','Color','w'); imshow(rgb_clean); hold on;
    for kk=1:n70
        idROI = idxClean(kk);
        m = (L_roi==idROI);
        B = bwboundaries(m,'noholes');
        for b=1:numel(B)
            plot(B{b}(:,2), B{b}(:,1), 'c', 'LineWidth', 0.9);
        end
    end
    title(sprintf('Selected clean ROIs (n=%d) — %s', n70, sampleName), 'Interpreter','none');
    hold off;
    savefig(h, fullfile(outSpectra70Dir, sprintf('STEP4_SelectedROIs70_%s.fig', sampleName)));
    exportgraphics(h, fullfile(outSpectra70Dir, sprintf('STEP4_SelectedROIs70_%s.png', sampleName)), 'Resolution', 300);
    close(h);

    h = figure('Visible','off','Color','w');
    plot(wl, Spectra70.','LineWidth',0.8);
    xlabel('Wavelength (nm)'); ylabel('Reflectance');
    title(sprintf('Mean spectra (n=%d) — %s', n70, sampleName), 'Interpreter','none');
    grid on; box on;
    savefig(h, fullfile(outSpectra70Dir, sprintf('STEP4_MeanSpectra70_%s.fig', sampleName)));
    exportgraphics(h, fullfile(outSpectra70Dir, sprintf('STEP4_MeanSpectra70_%s.png', sampleName)), 'Resolution', 300);
    close(h);

    fprintf('  OK: %s\n', sampleName);
end

% ==========================================================
% STEP 5: Compile 70-ROI tables + per-sample means
% ==========================================================
fprintf('\nSTEP 5: compiling 70-ROI tables + per-sample means...\n');

files70 = dir(fullfile(outSpectra70Dir, 'Spectra70_CentralROIs_*.xlsx'));
if isempty(files70)
    warning('No Spectra70_CentralROIs_*.xlsx found. Nothing to compile.');
    return;
end

Tall = table();
Tmeans = table();

metaVars = {'Sensor','Variety','SamplingDate','SampleIndex','Treatment','Block','IrrigationRegime','VineyardFloorManagement','SampleName'};

for iF=1:numel(files70)
    fpath = fullfile(files70(iF).folder, files70(iF).name);
    T = readtable(fpath, 'Sheet','Spectra70', 'TextType','string');

    if isempty(Tall)
        Tall = T;
    else
        % align columns if needed
        vAll = Tall.Properties.VariableNames;
        if ~isequal(vAll, T.Properties.VariableNames)
            T = T(:, vAll);
        end
        Tall = [Tall; T]; %#ok<AGROW>
    end

    % mean per sample
    Tm = table();
    for mv=1:numel(metaVars)
        if ismember(metaVars{mv}, T.Properties.VariableNames)
            Tm.(metaVars{mv}) = T.(metaVars{mv})(1);
        end
    end
    Tm.N_ROIs = height(T);

    isSpec = startsWith(T.Properties.VariableNames,'nm_');
    X = double(T{:,isSpec});
    mu = mean(X,1,'omitnan');
    specNames = T.Properties.VariableNames(isSpec);
    for j=1:numel(specNames)
        Tm.(specNames{j}) = mu(j);
    end

    if isempty(Tmeans), Tmeans = Tm; else, Tmeans = [Tmeans; Tm]; end %#ok<AGROW>
end

writetable(Tall,   fullfile(outSpectra70Dir,'Spectra70_ALL_Samples_compiled.xlsx'), 'Sheet','All_70ROIs','WriteMode','overwritesheet');
writetable(Tmeans, fullfile(outSpectra70Dir,'Spectra70_MEANS_bySample.xlsx'),        'Sheet','Means',   'WriteMode','overwritesheet');

fprintf('OK: compilation completed.\n');
end

%% ====================== helpers ======================
function ensureDir(p)
if ~exist(p,'dir'), mkdir(p); end
end

function design = makeDesignTable()
design = table( ...
    (1:24).', ...
    ["RLA";"RL";"RT";"RA";"SA";"SL";"RT";"RA";"RL";"RLA";"SL";"SA";"RA";"RL";"RLA";"RT";"SL";"SA";"RLA";"RT";"RL";"RA";"SA";"SL"], ...
    ["I";"I";"I";"I";"I";"I";"II";"II";"II";"II";"II";"II";"III";"III";"III";"III";"III";"III";"IV";"IV";"IV";"IV";"IV";"IV"], ...
    ["Irrigated";"Irrigated";"Irrigated";"Irrigated";"Rainfed";"Rainfed";"Irrigated";"Irrigated";"Irrigated";"Irrigated";"Rainfed";"Rainfed"; ...
     "Irrigated";"Irrigated";"Irrigated";"Irrigated";"Rainfed";"Rainfed";"Irrigated";"Irrigated";"Irrigated";"Irrigated";"Rainfed";"Rainfed"], ...
    ["Open vegetation";"Tillage";"Double cover crop";"Alternate cover crop";"Alternate cover crop";"Tillage"; ...
     "Double cover crop";"Alternate cover crop";"Tillage";"Open vegetation";"Tillage";"Alternate cover crop"; ...
     "Alternate cover crop";"Tillage";"Open vegetation";"Double cover crop";"Tillage";"Alternate cover crop"; ...
     "Open vegetation";"Double cover crop";"Tillage";"Alternate cover crop";"Alternate cover crop";"Tillage"], ...
    'VariableNames', {'SampleIndex','Treatment','Block','IrrigationRegime','VineyardFloorManagement'} );
end

function [cube_cal, nb] = calibrateReflectance(cube, white, dark)
nrow = min([size(cube,1), size(white,1), size(dark,1)]);
ncol = min([size(cube,2), size(white,2), size(dark,2)]);
nb   = min([size(cube,3), size(white,3), size(dark,3)]);
cube  = cube(1:nrow,1:ncol,1:nb);
white = white(1:nrow,1:ncol,1:nb);
dark  = dark(1:nrow,1:ncol,1:nb);

cube_cal = (double(cube) - double(dark)) ./ max(eps, (double(white) - double(dark)));
cube_cal(~isfinite(cube_cal)) = 0;
cube_cal = max(0, min(cube_cal, 1));
end

function wl = getWavelengthsFromInfoOrHdr(info, hdrPath, nb)
wl = [];
try
    if isstruct(info)
        if isfield(info,'wavelength') && ~isempty(info.wavelength), wl = double(info.wavelength(:)); end
        if isempty(wl) && isfield(info,'wavelengths') && ~isempty(info.wavelengths), wl = double(info.wavelengths(:)); end
    end
catch
    wl = [];
end

if (isempty(wl) || numel(wl)<nb) && exist(hdrPath,'file')
    try
        txt = fileread(hdrPath);
        tok = regexp(txt,'(?is)\bwavelength\s*=\s*\{(.*?)\}','tokens','once');
        if ~isempty(tok)
            nums = regexp(tok{1}, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
            wl = str2double(nums(:)); wl = wl(isfinite(wl));
        end
    catch
        wl = [];
    end
end

if isempty(wl)
    wl = linspace(950, 1720, nb).';
else
    wl = wl(:);
    if numel(wl) >= nb, wl = wl(1:nb);
    else
        wl = [wl; linspace(wl(end), wl(end)+10, nb-numel(wl)).'];
    end
end
end

function wl_lbls = makeNmLabels_dec2(wl)
w = double(wl(:));
lbl = arrayfun(@(x) sprintf('nm_%.2f',x), w, 'UniformOutput', false);
lbl = strrep(lbl,'.','_');
lbl = matlab.lang.makeValidName(lbl);
wl_lbls = matlab.lang.makeUniqueStrings(lbl(:));
end

function [Sensor, Variety, DateStr, SampleIndex, Treatment, Block, IrrigationRegime, VineyardFloorManagement] = parseMetaFromSampleName(sampleName, design)
tokens = split(string(sampleName), "_");

Sensor = ""; if numel(tokens)>=1, Sensor = tokens(1); end

Variety = "";
if numel(tokens)>=2
    vUp = upper(strtrim(tokens(2)));
    if vUp=="GN", Variety="Garnacha Tinta";
    elseif vUp=="GB", Variety="Garnacha Blanca";
    else, Variety=tokens(2);
    end
end

DateStr = "";
if numel(tokens)>=3
    tok = tokens(3);
    if strlength(tok)==4 && all(isstrprop(char(tok),'digit'))
        dd = str2double(extractBetween(tok,1,2));
        mm = str2double(extractBetween(tok,3,4));
        if isfinite(dd)&&isfinite(mm), DateStr = string(sprintf('%02d/%02d',dd,mm));
        else, DateStr = tok;
        end
    else
        DateStr = tok;
    end
end

SampleIndex = NaN;
if numel(tokens)>=1, SampleIndex = str2double(tokens(end)); end

Treatment=""; Block=""; IrrigationRegime=""; VineyardFloorManagement="";
if isfinite(SampleIndex)
    hit = design(design.SampleIndex==SampleIndex,:);
    if ~isempty(hit)
        Treatment = hit.Treatment(1);
        Block = hit.Block(1);
        IrrigationRegime = hit.IrrigationRegime(1);
        VineyardFloorManagement = hit.VineyardFloorManagement(1);
    end
end
end
