function fx10_batch_berry_extraction_surgical_precision()
% fx10_batch_berry_extraction_surgical_precision
%
% FX10 batch pipeline for berry-level spectral extraction.
%
% Strategy:
% 1. SEGMENTATION
%    - SAM on the original image (maximum detection).
%    - Strict NDVI live mask, without hole filling.
%
% 2. PURITY CUT
%    - Each SAM-detected object is intersected with the live mask.
%    - Background pixels incorrectly included by SAM are physically removed.
%
% 3. CLEANING AND SELECTION
%    - Solidity filter > 0.88.
%    - Size filter (10th percentile).
%    - Final 1-pixel erosion before spectral extraction (edge avoidance).
%    - Statistical selection 90 -> 70.
%
% Requirements:
%   - imsegsam
%   - Image Processing Toolbox

close all; clc;

%% --- WORKING ROOT ---
ruta_principal = [pwd filesep];
fprintf('Searching for samples and helper functions in: %s\n', ruta_principal);

%% --- SAM availability ---
if ~exist('imsegsam','file')
    error(['imsegsam was not found.\n' ...
           'Install "Image Processing Toolbox Model for Segment Anything Model" (Add-Ons).']);
end

%% --- Bind ENVI readers ---
read_info = @read_envi_info;
read_cube = @read_envi_cube;
if exist('enviinfo_local','file'), read_info = @enviinfo_local; end
if exist('enviread_local','file'), read_cube = @enviread_local; end

%% --- OUTPUT FOLDERS ---
carpeta_salida12 = fullfile(ruta_principal, 'resultados_PASO_1_2_v34_SURGICAL');
if ~exist(carpeta_salida12, 'dir'); mkdir(carpeta_salida12); end
carpeta_90 = fullfile(ruta_principal, 'resultados_PASO_3_ROIs_full_90');
carpeta_70 = fullfile(ruta_principal, 'resultados_PASO_3_ROIs_full_70');
if ~exist(carpeta_90,'dir'); mkdir(carpeta_90); end
if ~exist(carpeta_70,'dir'); mkdir(carpeta_70); end

%% --- AUTOMATIC SAMPLE DISCOVERY ---
d = dir(ruta_principal);
es_muestra = @(x) x.isdir && startsWith(x.name, 'FX10_') && ~contains(lower(x.name), 'resultado');
carpetas_muestra = d(arrayfun(es_muestra, d));
if isempty(carpetas_muestra)
    error('No sample folders (FX10_...) were found in: %s', ruta_principal);
end
fprintf('Found %d FX10 samples:\n', numel(carpetas_muestra));

%% ---------- PARAMETERS ----------
% 1. SAM
SAM_GRID_SIZE = [32 32];
SAM_SCORE     = 0.40;

% 2. NDVI (live mask)
WL_NIR = 800;
WL_RED = 670;
STRICT_NDVI_THRESHOLD = 0.38;

% 3. Quality filters
MIN_LABEL_AREA   = 150;   % Absolute minimum area (before percentile filter)
MIN_SOLIDITY     = 0.88;  % Compactness (avoids irregular shapes)

% 4. Statistical size filter
PERCENTILE_CUT = 10; % Remove the smallest 10%

% 5. Spectral selection
N_KEEP_90 = 90;
N_KEEP_70 = 70;

%% --- EXPERIMENTAL DESIGN ---
design = table( ...
    (1:24).', ...
    ["RLA";"RL";"RT";"RA";"SA";"SL";"RT";"RA";"RL";"RLA";"SL";"SA";"RA";"RL";"RLA";"RT";"SL";"SA";"RLA";"RT";"RL";"RA";"SA";"SL"], ...
    ["I";"I";"I";"I";"I";"I";"II";"II";"II";"II";"II";"II";"III";"III";"III";"III";"III";"III";"IV";"IV";"IV";"IV";"IV";"IV"], ...
    ["Riego";"Riego";"Riego";"Riego";"Secano";"Secano";"Riego";"Riego";"Riego";"Riego";"Secano";"Secano";"Riego";"Riego";"Riego";"Riego";"Secano";"Secano";"Riego";"Riego";"Riego";"Riego";"Secano";"Secano"], ...
    ["Vegetación abierta";"Laboreo";"Cubierta doble";"Cubierta alterna";"Cubierta alterna";"Laboreo"; ...
     "Cubierta doble";"Cubierta alterna";"Laboreo";"Vegetación abierta";"Laboreo";"Cubierta alterna"; ...
     "Cubierta alterna";"Laboreo";"Vegetación abierta";"Cubierta doble";"Laboreo";"Cubierta alterna"; ...
     "Vegetación abierta";"Cubierta doble";"Laboreo";"Cubierta alterna";"Cubierta alterna";"Laboreo"], ...
    'VariableNames', {'SampleIndex','Tratamiento','Bloque','Riego','Manejo_cultivo'} );

%% ========= MAIN LOOP =========
for k = 1:numel(carpetas_muestra)
    nombre_muestra = carpetas_muestra(k).name;
    ruta_captura   = fullfile(ruta_principal, nombre_muestra, 'capture');
    fprintf('\n=====================================================\n');
    fprintf('>> Processing sample: %s (%d/%d)\n', nombre_muestra, k, numel(carpetas_muestra));

    %% --- Loading and calibration ---
    try
        hdr_muestra = fullfile(ruta_captura, [nombre_muestra, '.hdr']);
        raw_muestra = fullfile(ruta_captura, [nombre_muestra, '.raw']);
        hdr_white   = fullfile(ruta_captura, ['WHITEREF_', nombre_muestra, '.hdr']);
        raw_white   = fullfile(ruta_captura, ['WHITEREF_', nombre_muestra, '.raw']);
        hdr_dark    = fullfile(ruta_captura, ['DARKREF_',  nombre_muestra, '.hdr']);
        raw_dark    = fullfile(ruta_captura, ['DARKREF_',  nombre_muestra, '.raw']);

        if ~exist(hdr_white,'file'), error('WHITE reference not found'); end

        info  = read_info(hdr_muestra);
        cube  = read_cube(info,  raw_muestra);
        infoW = read_info(hdr_white);   white = read_cube(infoW, raw_white);
        infoD = read_info(hdr_dark);    dark  = read_cube(infoD, raw_dark);

        nrow = min([size(cube,1), size(white,1), size(dark,1)]);
        ncol = min([size(cube,2), size(white,2), size(dark,2)]);
        nb   = min([size(cube,3), size(white,3), size(dark,3)]);
        cube  = cube (1:nrow, 1:ncol, 1:nb);
        white = white(1:nrow, 1:ncol, 1:nb);
        dark  = dark (1:nrow, 1:ncol, 1:nb);
        cube_cal = (double(cube) - double(dark)) ./ max(eps, (double(white) - double(dark)));
        cube_cal(~isfinite(cube_cal)) = 0;
        cube_cal = max(0, min(cube_cal, 1));
    catch ME
        fprintf('  ERROR loading data: %s\n', ME.message);
        continue;
    end

    wl = getWavelengths_FX_hdrFirst(info, hdr_muestra, nb);

    %% --- Visible image (RGB) ---
    png_root = fullfile(ruta_principal, nombre_muestra, [nombre_muestra '.png']);
    if exist(png_root,'file')
        I_vis = imread(png_root);
        srcVis = "PNG_sample_root";
    else
        [idxR, idxG, idxB] = pickRGBbands(wl, nb);
        I_vis = cat(3, cube_cal(:,:,idxR), cube_cal(:,:,idxG), cube_cal(:,:,idxB));
        I_vis = im2uint8(imadjust(I_vis, stretchlim(I_vis,[0.01 0.99])));
        srcVis = "pseudoRGB_from_cube";
    end
    rows = min(size(I_vis,1), nrow);
    cols = min(size(I_vis,2), ncol);
    I_vis    = I_vis(1:rows, 1:cols, :);
    cube_cal = cube_cal(1:rows, 1:cols, :);

    %% =======================================================
    %  STEP 1: LIVE MASK (strict NDVI, no hole filling)
    %% =======================================================
    fprintf('\n>> STEP 1: Computing live mask (NDVI)...\n');
    idxNIR = nearestBand(wl, WL_NIR, nb);
    idxRED = nearestBand(wl, WL_RED, nb);
    RNIR = cube_cal(:,:,idxNIR);
    RRED = cube_cal(:,:,idxRED);
    NDVI = (RNIR - RRED) ./ max(eps, (RNIR + RRED));
    NDVI(~isfinite(NDVI)) = -1;

    % Strict live mask
    % No imfill is applied here. If there is a hole, it is considered
    % background or non-living shaded tissue.
    mask_live = NDVI > STRICT_NDVI_THRESHOLD;
    mask_live = bwareaopen(mask_live, 30);

    h1 = figure('Name','STEP 1: Live Mask','WindowState','maximized','Color','w');
    imshow(mask_live);
    title(sprintf('Live Mask (Pure NDVI) - %s', nombre_muestra), 'Interpreter','none');
    drawnow;
    exportgraphics(h1, fullfile(carpeta_salida12, sprintf('PASO1_Mask_Live_%s.png', nombre_muestra)), 'Resolution', 250);

    %% =======================================================
    %  STEP 2: SAM + SURGICAL CUT
    %% =======================================================
    fprintf('\n>> STEP 2: SAM + cut + solidity filtering...\n');

    % Run SAM on the original image
    masksSAM = imsegsam(im2uint8(I_vis), ...
        'PointGridSize', SAM_GRID_SIZE, ...
        'ScoreThreshold', SAM_SCORE);
    Lsam_raw = labelmatrix(masksSAM);

    stats = regionprops(Lsam_raw, 'PixelIdxList');
    nlab = numel(stats);

    Lberries_temp = zeros(size(Lsam_raw),'uint16');
    new_id = 0;

    % --- SURGICAL LOOP ---
    for i = 1:nlab
        % 1. Get the original SAM object mask
        mask_obj = (Lsam_raw == i);

        % 2. CUT: intersect with live mask
        % Background pixels incorrectly included by SAM are removed.
        mask_cut = mask_obj & mask_live;

        if nnz(mask_cut) == 0, continue; end

        % If the cut splits the berry in two parts, keep the largest one
        mask_cut = keepLargestComponentDeterministic(mask_cut, mask_obj);

        % 3. Compute properties on the CUT object
        prop = regionprops(mask_cut, 'Area','Solidity');
        A = prop.Area;
        Sol = prop.Solidity;

        % 4. Quality filters
        if A < MIN_LABEL_AREA, continue; end
        if Sol < MIN_SOLIDITY, continue; end

        % Keep the cut object if it passes all filters
        new_id = new_id + 1;
        Lberries_temp(mask_cut) = new_id;
    end

    count_temp = new_id;
    fprintf('  Berries (pre-size filter): %d\n', count_temp);

    %% =======================================================
    %  STEP 2b: STATISTICAL SIZE FILTER (10th percentile)
    %% =======================================================
    if count_temp > 0
        statsT = regionprops(Lberries_temp, 'Area','PixelIdxList');
        allAreas = [statsT.Area];

        % 10th percentile threshold
        area_threshold_10 = prctile(allAreas, PERCENTILE_CUT);
        % Optional upper threshold (avoid oversized blobs)
        area_threshold_99 = prctile(allAreas, 99.5);

        fprintf('  -> Filtering: Area < %.0f px or > %.0f px\n', area_threshold_10, area_threshold_99);

        Lberries = zeros(size(Lberries_temp), 'uint16');
        new_id_final = 0;

        for k = 1:numel(statsT)
            A = statsT(k).Area;
            if A >= area_threshold_10 && A <= area_threshold_99
                new_id_final = new_id_final + 1;
                Lberries(statsT(k).PixelIdxList) = new_id_final;
            end
        end
        num_bayas = new_id_final;
    else
        Lberries = Lberries_temp;
        num_bayas = 0;
    end

    fprintf('  Final purified berries: %d\n', num_bayas);

    % --- FINAL FIGURES ---
    h2 = figure('Name','STEP 2: Final Berries','WindowState','maximized','Color','w');
    imshow(label2rgb(Lberries,'jet','w','shuffle'));
    title(sprintf('STEP 2: %d berries (Surgical Cut) - %s', num_bayas, nombre_muestra), 'Interpreter','none');
    drawnow;
    exportgraphics(h2, fullfile(carpeta_salida12, sprintf('PASO2_Bayas_Surgical_%s.png', nombre_muestra)), 'Resolution', 250);

    h2b = figure('Name','STEP 2: Contours','WindowState','maximized','Color','w');
    imshow(I_vis); hold on;
    for idb = 1:num_bayas
        m = (Lberries == idb);
        B = bwboundaries(m,'noholes');
        for kB = 1:numel(B)
            plot(B{kB}(:,2), B{kB}(:,1), 'c', 'LineWidth', 0.8);
        end
    end
    title(sprintf('STEP 2: Contours - %s', nombre_muestra), 'Interpreter','none');
    hold off; drawnow;
    exportgraphics(h2b, fullfile(carpeta_salida12, sprintf('PASO2_Contornos_Surgical_%s.png', nombre_muestra)), 'Resolution', 250);

    save(fullfile(carpeta_salida12, sprintf('labels_SAM_%s.mat', nombre_muestra)), 'Lberries');

    %% ==========================================================
    %  STEP 3: SPECTRAL EXTRACTION WITH EROSION
    %% ==========================================================
    fprintf('\n>> STEP 3: Spectral extraction (light erosion)... \n');

    Spectra_ALL = zeros(num_bayas, nb);

    for ii = 1:num_bayas
        mask_baya = (Lberries == ii);

        % Light erosion to remove the outermost edge
        mask_baya = imerode(mask_baya, strel('disk', 1));

        if nnz(mask_baya) == 0, continue; end

        pix = cube_cal(repmat(mask_baya, [1 1 nb]));
        pix = reshape(pix, [], nb);
        Spectra_ALL(ii,:) = mean(pix, 1);
    end

    %% ==========================================================
    %  STEP 4: HOMOGENEITY-BASED SELECTION (90 -> 70)
    %% ==========================================================
    fprintf('\n>> STEP 4: Statistical selection (90 -> 70)...\n');

    % --- PHASE A: 90 vs global median ---
    Mediana_Global = median(Spectra_ALL, 1);
    Dist_Global    = sqrt(sum((Spectra_ALL - Mediana_Global).^2, 2));
    [~, sort_idx_global] = sort(Dist_Global, 'ascend');

    n_90 = min(num_bayas, N_KEEP_90);
    idx_90 = sort_idx_global(1:n_90);
    Spectra_90 = Spectra_ALL(idx_90, :);
    Dist_90_fromGlobal = Dist_Global(idx_90);

    % --- PHASE B: 70 vs median of the selected 90 ---
    Mediana_90 = median(Spectra_90, 1);
    Dist_Local = sqrt(sum((Spectra_90 - Mediana_90).^2, 2));
    [~, sort_idx_local] = sort(Dist_Local, 'ascend');

    n_70 = min(n_90, N_KEEP_70);
    idx_70_local = sort_idx_local(1:n_70);
    idx_70_final = idx_90(idx_70_local);
    Spectra_70   = Spectra_90(idx_70_local, :);
    Dist_70_fromLocal = Dist_Local(idx_70_local);

    % --- EXPORT ---
    wl_lbls = makeNmLabels_dec2(wl);

    % EXCEL 90
    SampleName = repmat(string(nombre_muestra), n_90, 1);
    Tspec90 = table(SampleName, idx_90, Dist_90_fromGlobal, 'VariableNames', {'SampleName','ROI_ID_Original','Dist_to_GlobalMedian'});
    for j = 1:nb, Tspec90.(wl_lbls{j}) = Spectra_90(:,j); end
    writetable(Tspec90, fullfile(carpeta_90, sprintf('Spectra90_FullROI_%s.xlsx', nombre_muestra)), 'Sheet','Spectra90','WriteMode','overwritesheet');

    % EXCEL 70
    [Camara, Variedad, Fecha, SampleIndex_val, Tratamiento_val, Bloque_val, Riego_val, Manejo_val] = ...
        parseMetaFromSampleName(nombre_muestra, design);

    SampleName70 = repmat(string(nombre_muestra), n_70, 1);
    Tspec70 = table( ...
        repmat(string(Camara), n_70, 1), ...
        repmat(string(Variedad), n_70, 1), ...
        repmat(string(Fecha), n_70, 1), ...
        repmat(double(SampleIndex_val), n_70, 1), ...
        repmat(string(Tratamiento_val), n_70, 1), ...
        repmat(string(Bloque_val), n_70, 1), ...
        repmat(string(Riego_val), n_70, 1), ...
        repmat(string(Manejo_val), n_70, 1), ...
        SampleName70, idx_70_final, Dist_70_fromLocal, ...
        'VariableNames', {'Camara','Variedad','Fecha','SampleIndex','Tratamiento','Bloque','Riego','Manejo_cultivo','SampleName','ROI_ID_Original','Dist_to_Median90'} );
    for j = 1:nb, Tspec70.(wl_lbls{j}) = Spectra_70(:,j); end
    writetable(Tspec70, fullfile(carpeta_70, sprintf('Spectra70_FullROI_%s.xlsx', nombre_muestra)), 'Sheet','Spectra70','WriteMode','overwritesheet');

    % --- QC FIGURES ---
    mask_70_vis = ismember(Lberries, idx_70_final);
    I_nobg_70 = applyWhiteBackground(I_vis, mask_70_vis);

    hfig70 = figure('Name','STEP 4: Selected 70 ROIs','WindowState','maximized','Color','w');
    imshow(I_nobg_70); hold on;
    for k_idx = 1:n_70
        idr = idx_70_final(k_idx);
        m = (Lberries == idr);
        B = bwboundaries(m,'noholes');
        for kB = 1:numel(B), plot(B{kB}(:,2), B{kB}(:,1), 'c', 'LineWidth', 1.2); end
    end
    title(sprintf('Top 70 Clean Berries (Surgical Cut) - %s', nombre_muestra), 'Interpreter','none');
    hold off; drawnow;
    file_fig70 = fullfile(carpeta_70, sprintf('PASO3_ROIsSelected70_%s.fig', nombre_muestra));
    file_jpg70 = fullfile(carpeta_70, sprintf('PASO3_ROIsSelected70_%s.jpg', nombre_muestra));
    savefig(hfig70, file_fig70);
    exportgraphics(hfig70, file_jpg70, 'Resolution', 300);

    hspec70 = figure('Name','STEP 4: Spectra 70','WindowState','maximized','Color','w');
    plot(wl, Spectra_70.', 'LineWidth', 0.8);
    xlabel('nm'); ylabel('Reflectance');
    title(sprintf('Mean spectra 70 (n=%d) - %s', n_70, nombre_muestra), 'Interpreter','none');
    grid on; box on; drawnow;
    file_figS70 = fullfile(carpeta_70, sprintf('PASO3_MeanSpectra70_%s.fig', nombre_muestra));
    file_jpgS70 = fullfile(carpeta_70, sprintf('PASO3_MeanSpectra70_%s.jpg', nombre_muestra));
    savefig(hspec70, file_figS70);
    exportgraphics(hspec70, file_jpgS70, 'Resolution', 300);

    close all;
    fprintf('  [OK] %s\n', nombre_muestra);
end
fprintf('\n>> Per-sample processing completed.\n');

%% ==========================================================
%  STEP 5: COMPILE
%% ==========================================================
fprintf('\n>> STEP 5: Compiling...\n');
files70 = dir(fullfile(carpeta_70, 'Spectra70_FullROI_*.xlsx'));
if isempty(files70)
    warning('No Excel files were found.'); return;
end
Tall = table(); Tmeans = table();
for iF = 1:numel(files70)
    fpath = fullfile(files70(iF).folder, files70(iF).name);
    try, T = readtable(fpath, 'Sheet', 'Spectra70', 'TextType','string'); catch, T = readtable(fpath, 'Sheet', 'Spectra70'); end
    if isempty(Tall), Tall = T; else
        if ~isequal(Tall.Properties.VariableNames, T.Properties.VariableNames)
             T = T(:, Tall.Properties.VariableNames);
        end
        Tall = [Tall; T]; %#ok<AGROW>
    end

    Tm = table();
    metaVars = {'Camara','Variedad','Fecha','SampleIndex','Tratamiento','Bloque','Riego','Manejo_cultivo','SampleName'};
    for mv = 1:numel(metaVars), Tm.(metaVars{mv}) = T.(metaVars{mv})(1); end
    Tm.N_ROIs = height(T);
    isSpec = startsWith(T.Properties.VariableNames, 'nm_');
    X = double(T{:, isSpec});
    mu = mean(X, 1, 'omitnan');
    sNames = T.Properties.VariableNames(isSpec);
    for j=1:numel(sNames), Tm.(sNames{j}) = mu(j); end

    if isempty(Tmeans), Tmeans = Tm; else, Tmeans = [Tmeans; Tm]; %#ok<AGROW>
    end
end
writetable(Tall, fullfile(carpeta_70, 'Spectra70_ALL_Samples_compiled.xlsx'), 'Sheet','All_70ROIs','WriteMode','overwritesheet');
writetable(Tmeans, fullfile(carpeta_70, 'Spectra70_MEANS_bySample.xlsx'), 'Sheet','Means','WriteMode','overwritesheet');
fprintf('  [OK] Compiled files created.\n');
end

%% ===================== HELPERS =====================
function wl = getWavelengths_FX_hdrFirst(info, hdrPath, nb)
wl = [];
try
    if isstruct(info)
        if isfield(info,'wavelength') && ~isempty(info.wavelength), wl = double(info.wavelength(:));
        elseif isfield(info,'wavelengths') && ~isempty(info.wavelengths), wl = double(info.wavelengths(:));
        elseif isfield(info,'Wavelength') && ~isempty(info.Wavelength), wl = double(info.Wavelength(:)); end
    end
catch, wl = []; end
if (isempty(wl) || numel(wl) < nb) && exist(hdrPath,'file')
    try
        txt = fileread(hdrPath);
        tok = regexp(txt, '(?is)\bwavelength\s*=\s*\{(.*?)\}', 'tokens', 'once');
        if ~isempty(tok)
            nums = regexp(tok{1}, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
            wl = str2double(nums(:)); wl = wl(isfinite(wl));
        end
    catch, wl = []; end
end
if isempty(wl) || numel(wl) < nb, if nb == 224, wl = fx10_wavelengths_224(); end; end
if isempty(wl), wl = linspace(400, 1000, nb).'; end
if numel(wl) >= nb, wl = wl(1:nb); else, wl = [wl(:); linspace(wl(end), 1000, nb-numel(wl)).']; end
end

function [idxR, idxG, idxB] = pickRGBbands(wl, nb)
idxR = nearestBand(wl, 650, nb); idxG = nearestBand(wl, 550, nb); idxB = nearestBand(wl, 450, nb);
end

function idx = nearestBand(wl, target, nb)
if isempty(wl) || numel(wl) < nb, idx = max(1, min(nb, round((target-400)/(1000-400)*(nb-1)+1))); return; end
[~, idx] = min(abs(wl - target)); idx = max(1, min(nb, idx));
end
function mask_out = keepLargestComponentDeterministic(mask_in, ref_mask)
% Keep the largest connected component.
% If multiple components share the same maximum area, break the tie
% deterministically using the centroid distance to the original SAM object.
%
% Inputs:
%   mask_in  - binary mask after NDVI-based trimming
%   ref_mask - original SAM object mask before trimming
%
% Output:
%   mask_out - binary mask containing one selected component only

mask_out = false(size(mask_in));

if ~any(mask_in(:))
    return;
end

cc = bwconncomp(mask_in);
if cc.NumObjects == 1
    mask_out = mask_in;
    return;
end

areas = cellfun(@numel, cc.PixelIdxList);
maxArea = max(areas);
cand = find(areas == maxArea);

% Unique largest component
if numel(cand) == 1
    mask_out(cc.PixelIdxList{cand}) = true;
    return;
end

% Deterministic tie-break: closest centroid to the original SAM object
if nargin < 2 || isempty(ref_mask) || ~any(ref_mask(:))
    ref_mask = mask_in;
end

refStats = regionprops(ref_mask, 'Centroid');
refCentroid = refStats(1).Centroid;

compStats = regionprops(cc, 'Centroid');
distances = inf(numel(cand), 1);

for ii = 1:numel(cand)
    c = compStats(cand(ii)).Centroid;
    distances(ii) = hypot(c(1) - refCentroid(1), c(2) - refCentroid(2));
end

minDist = min(distances);
best = cand(distances == minDist);

% Final deterministic fallback in the unlikely event of another tie
best = best(1);

mask_out(cc.PixelIdxList{best}) = true;
end

function Iout = applyWhiteBackground(Iin, maskKeep)
Iout = Iin; if isempty(maskKeep), return; end
if size(Iout,3) == 1, if isinteger(Iout), Iout(~maskKeep) = intmax(class(Iout)); else, Iout(~maskKeep) = 1; end
else, for c = 1:size(Iout,3), chan = Iout(:,:,c); if isinteger(chan), chan(~maskKeep) = intmax(class(chan)); else, chan(~maskKeep) = 1; end; Iout(:,:,c) = chan; end; end
end

function wl_lbls = makeNmLabels_dec2(wl)
w = double(wl(:)); lbl = arrayfun(@(x) sprintf('nm_%.2f', x), w, 'UniformOutput', false);
lbl = strrep(lbl, '.', '_'); lbl = matlab.lang.makeUniqueStrings(lbl); wl_lbls = lbl(:);
end

function [Camara, Variedad, Fecha, SampleIndex_val, Tratamiento_val, Bloque_val, Riego_val, Manejo_val] = parseMetaFromSampleName(nombre_muestra, design)
tokens = split(string(nombre_muestra), "_"); Camara = ""; if numel(tokens) >= 1, Camara = tokens(1); end
Variedad = ""; if numel(tokens) >= 2, vUp = upper(strtrim(tokens(2))); if vUp == "GN", Variedad = "Garnacha tinta"; elseif vUp == "GB", Variedad = "Garnacha blanca"; else, Variedad = tokens(2); end; end
Fecha = ""; if numel(tokens) >= 3, tokFecha = tokens(3); if strlength(tokFecha)==4 && all(isstrprop(char(tokFecha),'digit')), dd = str2double(extractBetween(tokFecha,1,2)); mm = str2double(extractBetween(tokFecha,3,4)); if isfinite(dd) && isfinite(mm), Fecha = string(sprintf('%02d/%02d', dd, mm)); else, Fecha = tokFecha; end; else, Fecha = tokFecha; end; end
SampleIndex_val = NaN; if numel(tokens) >= 1, SampleIndex_val = str2double(tokens(end)); end
Tratamiento_val=""; Bloque_val=""; Riego_val=""; Manejo_val="";
if isfinite(SampleIndex_val), hit = design(design.SampleIndex == SampleIndex_val, :); if ~isempty(hit), Tratamiento_val = hit.Tratamiento(1); Bloque_val = hit.Bloque(1); Riego_val = hit.Riego(1); Manejo_val = hit.Manejo_cultivo(1); end; end
end

function wl = fx10_wavelengths_224()
wl = [400.67; 403.35; 406.03; 408.71; 411.38; 414.06; 416.74; 419.42; 422.10; 424.78; 427.46; 430.13; 432.81; 435.49; 438.17; 440.85; 443.53; 446.21; 448.88; 451.56; 454.24; 456.92; 459.60; 462.28; 464.96; 467.63; 470.31; 472.99; 475.67; 478.35; 481.03; 483.71; 486.38; 489.06; 491.74; 494.42; 497.10; 499.78; 502.46; 505.13; 507.81; 510.49; 513.17; 515.85; 518.53; 521.21; 523.88; 526.56; 529.24; 531.92; 534.60; 537.28; 539.96; 542.63; 545.31; 547.99; 550.67; 553.35; 556.03; 558.71; 561.38; 564.06; 566.74; 569.42; 572.10; 574.78; 577.46; 580.13; 582.81; 585.49; 588.17; 590.85; 593.53; 596.21; 598.88; 601.56; 604.24; 606.92; 609.60; 612.28; 614.96; 617.63; 620.31; 622.99; 625.67; 628.35; 631.03; 633.71; 636.38; 639.06; 641.74; 644.42; 647.10; 649.78; 652.46; 655.13; 657.81; 660.49; 663.17; 665.85; 668.53; 671.21; 673.88; 676.56; 679.24; 681.92; 684.60; 687.28; 689.96; 692.63; 695.31; 697.99; 700.67; 703.35; 706.03; 708.71; 711.38; 714.06; 716.74; 719.42; 722.10; 724.78; 727.46; 730.13; 732.81; 735.49; 738.17; 740.85; 743.53; 746.21; 748.88; 751.56; 754.24; 756.92; 759.60; 762.28; 764.96; 767.63; 770.31; 772.99; 775.67; 778.35; 781.03; 783.71; 786.38; 789.06; 791.74; 794.42; 797.10; 799.78; 802.46; 805.13; 807.81; 810.49; 813.17; 815.85; 818.53; 821.21; 823.88; 826.56; 829.24; 831.92; 834.60; 837.28; 839.96; 842.63; 845.31; 847.99; 850.67; 853.35; 856.03; 858.71; 861.38; 864.06; 866.74; 869.42; 872.10; 874.78; 877.46; 880.13; 882.81; 885.49; 888.17; 890.85; 893.53; 896.21; 898.88; 901.56; 904.24; 906.92; 909.60; 912.28; 914.96; 917.63; 920.31; 922.99; 925.67; 928.35; 931.03; 933.71; 936.38; 939.06; 941.74; 944.42; 947.10; 949.78; 952.46; 955.13; 957.81; 960.49; 963.17; 965.85; 968.53; 971.21; 973.88; 976.56; 979.24; 981.92; 984.60; 987.28; 989.96; 992.63; 995.31; 997.99];
end