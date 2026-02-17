function Automatizacion_bayas_FX10_por_lote_v34_SURGICAL_PRECISION_compile()
% Automatizacion_bayas_FX10_por_lote_v34_SURGICAL_PRECISION_compile
%
% FX10 (VIS–NIR) batch pipeline for berry-level spectral extraction.
%
% Pipeline (algorithmic intent preserved):
%  1) Reflectance calibration using sample-specific WHITE/DARK references.
%  2) Strict NDVI live mask (no hole filling).
%  3) SAM segmentation on the visual image; "surgical cut" = SAM object ∩ NDVI mask.
%  4) Geometric filtering (min area + min solidity) + area trimming (P10–P99.5).
%  5) Mean spectra per berry ROI after 1-pixel erosion.
%  6) Spectral QC selection: keep 90 closest to the global median, then 70 closest to the 90-median.
%  7) Export per-sample Excel tables + QC figures, then compile across samples.
%
% Expected structure (run from experiment root):
%   ROOT/
%     FX10_.../
%       capture/
%         <sample>.hdr/.raw
%         WHITEREF_<sample>.hdr/.raw
%         DARKREF_<sample>.hdr/.raw
%       <sample>.png   (optional; if absent, pseudo-RGB is built from cube)
%
% Outputs (created in ROOT):
%   - results_STEP1_2_FX10_v34_SURGICAL/
%   - results_Spectra90_FX10_FullROI/
%   - results_Spectra70_FX10_FullROI/
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
outStep12Dir     = fullfile(rootDir, 'results_STEP1_2_FX10_v34_SURGICAL');
outSpectra90Dir  = fullfile(rootDir, 'results_Spectra90_FX10_FullROI');
outSpectra70Dir  = fullfile(rootDir, 'results_Spectra70_FX10_FullROI');
ensureDir(outStep12Dir);
ensureDir(outSpectra90Dir);
ensureDir(outSpectra70Dir);

% --- Discover samples ---
d = dir(rootDir);
isSample = @(x) x.isdir && startsWith(x.name,'FX10_') && ~contains(lower(x.name),'resultado');
sampleDirs = d(arrayfun(isSample,d));
if isempty(sampleDirs)
    error('No FX10_* sample folders found in: %s', rootDir);
end
fprintf('FX10 samples found: %d\n', numel(sampleDirs));

% ---------- Parameters ----------
SAM_GRID_SIZE = [32 32];
SAM_SCORE     = 0.40;

WL_NIR = 800;
WL_RED = 670;
NDVI_THRESHOLD_STRICT = 0.38;

MIN_LABEL_AREA = 150;
MIN_SOLIDITY   = 0.88;

PERCENTILE_CUT_LOW  = 10;
PERCENTILE_CUT_HIGH = 99.5;

N_KEEP_90 = 90;
N_KEEP_70 = 70;

% --- Experimental design mapping (SampleIndex -> factors) ---
design = makeDesignTable();

% ========= Main loop =========
for k = 1:numel(sampleDirs)
    sampleName = sampleDirs(k).name;
    captureDir = fullfile(rootDir, sampleName, 'capture');

    fprintf('\n=====================================================\n');
    fprintf('Processing FX10 sample: %s (%d/%d)\n', sampleName, k, numel(sampleDirs));

    % --- Paths ---
    hdr_sample = fullfile(captureDir, [sampleName '.hdr']);
    raw_sample = fullfile(captureDir, [sampleName '.raw']);
    hdr_white  = fullfile(captureDir, ['WHITEREF_' sampleName '.hdr']);
    raw_white  = fullfile(captureDir, ['WHITEREF_' sampleName '.raw']);
    hdr_dark   = fullfile(captureDir, ['DARKREF_'  sampleName '.hdr']);
    raw_dark   = fullfile(captureDir, ['DARKREF_'  sampleName '.raw']);

    if ~exist(hdr_white,'file')
        fprintf('  ERROR: WHITE reference not found for %s\n', sampleName);
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

    % --- Visual image (PNG preferred) ---
    pngPath = fullfile(rootDir, sampleName, [sampleName '.png']);
    if exist(pngPath,'file')
        I_vis = imread(pngPath);
        visSource = "PNG";
    else
        [idxR, idxG, idxB] = pickRGBbands(wl, nb);
        I_vis = cat(3, cube_cal(:,:,idxR), cube_cal(:,:,idxG), cube_cal(:,:,idxB));
        I_vis = im2uint8(imadjust(I_vis, stretchlim(I_vis,[0.01 0.99])));
        visSource = "pseudoRGB";
    end

    % Align sizes (defensive)
    rows = min(size(I_vis,1), size(cube_cal,1));
    cols = min(size(I_vis,2), size(cube_cal,2));
    I_vis    = I_vis(1:rows,1:cols,:);
    cube_cal = cube_cal(1:rows,1:cols,:);

    % =======================================================
    % STEP 1: Strict NDVI live mask
    % =======================================================
    fprintf('  STEP 1: NDVI strict live mask (%s)\n', visSource);
    idxNIR = nearestBand(wl, WL_NIR, nb);
    idxRED = nearestBand(wl, WL_RED, nb);

    RNIR = cube_cal(:,:,idxNIR);
    RRED = cube_cal(:,:,idxRED);
    NDVI = (RNIR - RRED) ./ max(eps, (RNIR + RRED));
    NDVI(~isfinite(NDVI)) = -1;

    mask_live = NDVI > NDVI_THRESHOLD_STRICT;
    mask_live = bwareaopen(mask_live, 30);

    h = figure('Visible','off','Color','w'); imshow(mask_live);
    title(sprintf('STEP 1: NDVI live mask — %s', sampleName), 'Interpreter','none');
    exportgraphics(h, fullfile(outStep12Dir, sprintf('STEP1_NDVI_LiveMask_%s.png', sampleName)), 'Resolution', 250);
    close(h);

    % =======================================================
    % STEP 2: SAM segmentation + surgical cut + geometry filters
    % =======================================================
    fprintf('  STEP 2: SAM segmentation + surgical cut\n');

    masksSAM = imsegsam(im2uint8(I_vis), ...
        'PointGridSize', SAM_GRID_SIZE, ...
        'ScoreThreshold', SAM_SCORE);
    Lsam_raw = labelmatrix(masksSAM);

    stats = regionprops(Lsam_raw,'PixelIdxList');
    nlab = numel(stats);

    Lberries_temp = zeros(size(Lsam_raw),'uint16');
    new_id = 0;

    for i = 1:nlab
        objMask = (Lsam_raw == i);
        cutMask = objMask & mask_live;
        if nnz(cutMask)==0, continue; end

        cutMask = bwareafilt(cutMask, 1);
        p = regionprops(cutMask,'Area','Solidity');
        if isempty(p), continue; end

        if p.Area < MIN_LABEL_AREA, continue; end
        if p.Solidity < MIN_SOLIDITY, continue; end

        new_id = new_id + 1;
        Lberries_temp(cutMask) = new_id;
    end

    if new_id == 0
        fprintf('  WARNING: No berries detected after SAM+cut+filters for %s\n', sampleName);
        continue;
    end

    % =======================================================
    % STEP 2b: Area trimming (P10–P99.5)
    % =======================================================
    st = regionprops(Lberries_temp,'Area','PixelIdxList');
    areas = [st.Area];
    thr_low  = prctile(areas, PERCENTILE_CUT_LOW);
    thr_high = prctile(areas, PERCENTILE_CUT_HIGH);

    Lberries = zeros(size(Lberries_temp),'uint16');
    new_final = 0;
    for ii = 1:numel(st)
        A = st(ii).Area;
        if A >= thr_low && A <= thr_high
            new_final = new_final + 1;
            Lberries(st(ii).PixelIdxList) = new_final;
        end
    end
    numBerries = new_final;
    fprintf('  Final berries after trimming: %d\n', numBerries);

    % STEP2 figures
    h = figure('Visible','off','Color','w');
    imshow(label2rgb(Lberries,'jet','w','shuffle'));
    title(sprintf('STEP 2: %d berries — %s', numBerries, sampleName), 'Interpreter','none');
    exportgraphics(h, fullfile(outStep12Dir, sprintf('STEP2_Berries_Final_%s.png', sampleName)), 'Resolution', 250);
    close(h);

    h = figure('Visible','off','Color','w');
    imshow(I_vis); hold on;
    for idb = 1:numBerries
        m = (Lberries==idb);
        B = bwboundaries(m,'noholes');
        for b = 1:numel(B)
            plot(B{b}(:,2), B{b}(:,1), 'c', 'LineWidth', 0.8);
        end
    end
    title(sprintf('STEP 2: contours — %s', sampleName), 'Interpreter','none');
    hold off;
    exportgraphics(h, fullfile(outStep12Dir, sprintf('STEP2_Contours_%s.png', sampleName)), 'Resolution', 250);
    close(h);

    save(fullfile(outStep12Dir, sprintf('labels_SAM_%s.mat', sampleName)), 'Lberries');

    % =======================================================
    % STEP 3: Mean spectra per berry (1 px erosion)
    % =======================================================
    fprintf('  STEP 3: mean spectra per berry (edge-eroded)\n');

    Spectra_ALL = zeros(numBerries, nb);
    for ii = 1:numBerries
        m = (Lberries==ii);
        m = imerode(m, strel('disk',1));
        if nnz(m)==0, continue; end
        pix = cube_cal(repmat(m,[1 1 nb]));
        pix = reshape(pix,[],nb);
        Spectra_ALL(ii,:) = mean(pix,1);
    end

    % =======================================================
    % STEP 4: Selection 90 -> 70 by distance to median spectra
    % =======================================================
    fprintf('  STEP 4: spectral QC selection (90 -> 70)\n');

    globalMedian = median(Spectra_ALL,1);
    distGlobal   = sqrt(sum((Spectra_ALL - globalMedian).^2,2));
    [~, ordG]    = sort(distGlobal,'ascend');

    n90 = min(numBerries, N_KEEP_90);
    idx90 = ordG(1:n90);
    Spectra90 = Spectra_ALL(idx90,:);
    dist90_toGlobal = distGlobal(idx90);

    median90 = median(Spectra90,1);
    distLocal = sqrt(sum((Spectra90 - median90).^2,2));
    [~, ordL] = sort(distLocal,'ascend');

    n70 = min(n90, N_KEEP_70);
    idx70_local = ordL(1:n70);

    idx70_final = idx90(idx70_local);
    Spectra70   = Spectra90(idx70_local,:);
    dist70_toMedian90 = distLocal(idx70_local);

    % --- Export tables ---
    wlLabels = makeNmLabels_dec2(wl);

    % 90
    T90 = table(repmat(string(sampleName), n90,1), idx90(:), dist90_toGlobal(:), ...
        'VariableNames', {'SampleName','ROI_ID_Original','Dist_to_GlobalMedian'});
    for j=1:nb
        T90.(wlLabels{j}) = Spectra90(:,j);
    end
    writetable(T90, fullfile(outSpectra90Dir, sprintf('Spectra90_FullROI_%s.xlsx', sampleName)), ...
        'Sheet','Spectra90','WriteMode','overwritesheet');

    % 70 with metadata
    [Sensor, Variety, DateStr, SampleIndex, Treatment, Block, IrrigationRegime, VineyardFloorManagement] = ...
        parseMetaFromSampleName(sampleName, design);

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
        idx70_final(:), dist70_toMedian90(:), ...
        'VariableNames', {'Sensor','Variety','SamplingDate','SampleIndex','Treatment','Block','IrrigationRegime','VineyardFloorManagement', ...
                          'SampleName','ROI_ID_Original','Dist_to_Median90'} );
    for j=1:nb
        T70.(wlLabels{j}) = Spectra70(:,j);
    end
    writetable(T70, fullfile(outSpectra70Dir, sprintf('Spectra70_FullROI_%s.xlsx', sampleName)), ...
        'Sheet','Spectra70','WriteMode','overwritesheet');

    % --- QC figures ---
    mask70 = ismember(Lberries, idx70_final);
    I_whitebg = applyWhiteBackground(I_vis, mask70);

    h = figure('Visible','off','Color','w'); imshow(I_whitebg); hold on;
    for kk=1:n70
        idr = idx70_final(kk);
        m = (Lberries==idr);
        B = bwboundaries(m,'noholes');
        for b=1:numel(B)
            plot(B{b}(:,2), B{b}(:,1), 'c', 'LineWidth', 1.1);
        end
    end
    title(sprintf('Selected ROIs (n=%d) — %s', n70, sampleName), 'Interpreter','none');
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

files70 = dir(fullfile(outSpectra70Dir, 'Spectra70_FullROI_*.xlsx'));
if isempty(files70)
    warning('No Spectra70_FullROI_*.xlsx found. Nothing to compile.');
    return;
end

Tall = table();
Tmeans = table();

metaVars = {'Sensor','Variety','SamplingDate','SampleIndex','Treatment','Block','IrrigationRegime','VineyardFloorManagement','SampleName'};

for iF = 1:numel(files70)
    fpath = fullfile(files70(iF).folder, files70(iF).name);
    T = readtable(fpath, 'Sheet','Spectra70', 'TextType','string');

    % Append all rows
    if isempty(Tall)
        Tall = T;
    else
        if ~isequal(Tall.Properties.VariableNames, T.Properties.VariableNames)
            % align to Tall
            T = T(:, Tall.Properties.VariableNames);
        end
        Tall = [Tall; T]; %#ok<AGROW>
    end

    % Per-sample mean
    Tm = table();
    for mv = 1:numel(metaVars)
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
    wl = linspace(400,1000,nb).';
else
    wl = wl(:);
    if numel(wl) >= nb, wl = wl(1:nb);
    else
        wl = [wl; linspace(wl(end), 1000, nb-numel(wl)).'];
    end
end
end

function [idxR, idxG, idxB] = pickRGBbands(wl, nb)
idxR = nearestBand(wl, 650, nb);
idxG = nearestBand(wl, 550, nb);
idxB = nearestBand(wl, 450, nb);
end

function idx = nearestBand(wl, target, nb)
if isempty(wl) || numel(wl)<nb
    idx = max(1, min(nb, round((target-400)/(1000-400)*(nb-1)+1)));
    return;
end
[~,idx] = min(abs(wl-target));
idx = max(1, min(nb, idx));
end

function Iout = applyWhiteBackground(Iin, maskKeep)
Iout = Iin;
if isempty(maskKeep), return; end
if size(Iout,3)==1
    Iout(~maskKeep) = intmax(class(Iout));
else
    for c=1:size(Iout,3)
        ch = Iout(:,:,c);
        ch(~maskKeep) = intmax(class(ch));
        Iout(:,:,c) = ch;
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
