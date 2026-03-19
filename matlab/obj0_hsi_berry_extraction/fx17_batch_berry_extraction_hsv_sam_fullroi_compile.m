function fx17_batch_berry_extraction_hsv_sam_fullroi_compile()
% fx17_batch_berry_extraction_hsv_sam_fullroi_compile
%
% FX17 batch pipeline for berry-level spectral extraction.
%
% Based on the original v4_meta70 workflow, using full-berry ROIs with
% 2-pixel border erosion.
%
% Main features:
% (1) Variety mapping: "GN/GB" -> "Garnacha tinta / Garnacha blanca".
% (2) Final export:
%     - Compiled Excel file including all rows from all 70-ROI Excel files.
%     - Excel file with one mean spectrum per sample (one row per 70-ROI Excel).
%
% STEP 1 (HSV) and STEP 2 (SAM) are preserved.

close all; clc;

%% --- WORKING ROOT ---
ruta_principal = [pwd filesep];
fprintf('Searching for samples and helper functions in: %s\n', ruta_principal);

%% --- Check SAM availability ---
if ~exist('imsegsam','file')
    error(['imsegsam was not found.\n' ...
           'Install "Image Processing Toolbox Model for Segment Anything Model"\n' ...
           'from Add-Ons before running this script.']);
end

%% --- Bind ENVI readers ---
try
    read_info = @read_envi_info;
    read_cube = @read_envi_cube;
    if exist('enviinfo_local','file'), read_info = @enviinfo_local; end
    if exist('enviread_local','file'), read_cube = @enviread_local; end
catch ME
    fprintf('ERROR: Required helper functions were not found (read_envi_info.m, enviread_local.m, ...).\n');
    rethrow(ME);
end

%% --- OUTPUT FOLDERS (STEP 1–2) ---
carpeta_salida12 = fullfile(ruta_principal, 'resultados_PASO_1_SAM_AntiRaquis');
if ~exist(carpeta_salida12, 'dir'); mkdir(carpeta_salida12); end

%% --- AUTOMATIC SAMPLE DISCOVERY ---
d = dir(ruta_principal);
es_muestra = @(x) x.isdir && startsWith(x.name, 'FX17_') && ...
                         ~contains(lower(x.name), 'resultado');
carpetas_muestra = d(arrayfun(es_muestra, d));
if isempty(carpetas_muestra)
    error('No sample folders (FX17_...) were found in:%s\n%s', newline, ruta_principal);
end
fprintf('Found %d samples:\n', numel(carpetas_muestra));
for k_idx = 1:numel(carpetas_muestra)
    fprintf('  %s\n', carpetas_muestra(k_idx).name);
end

%% ---------- STEP 1 PARAMETERS (HSV) ----------
MIN_BERRY_AREA = 50;                  % small-noise removal
HUE_RANGE = [0.55 0.8];               % green grape range (anti-rachis)
MIN_SATURATION_THRESHOLD_HSV = 0.35;  % slightly stricter to remove paper/background
MAX_VALUE_THRESHOLD_HSV       = 0.85; % removes very bright regions

%% ---------- ROI PARAMETER ----------
ERODE_PIX = 2;   % ROI = full berry eroded by 2 px

%% --- EXPERIMENTAL DESIGN: SampleIndex -> factors (for metadata in Excel70) ---
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

%% ========= MAIN LOOP: PROCESS ALL FX17_ SAMPLES =========
for k = 1:numel(carpetas_muestra)

    nombre_muestra = carpetas_muestra(k).name;
    ruta_captura   = fullfile(ruta_principal, nombre_muestra, 'capture');

    fprintf('\n=====================================================\n');
    fprintf('>> Processing sample: %s (%d/%d)\n', ...
        nombre_muestra, k, numel(carpetas_muestra));
    fprintf('=====================================================\n');

    fprintf('\n>> STEP 1: HSV segmentation (anti-rachis)\n');
    fprintf('-----------------------------------------\n');

    %% --- Load and calibrate FX17 ---
    try
        hdr_muestra = fullfile(ruta_captura, [nombre_muestra, '.hdr']);
        raw_muestra = fullfile(ruta_captura, [nombre_muestra, '.raw']);

        hdr_white = fullfile(ruta_captura, ['WHITEREF_', nombre_muestra, '.hdr']);
        raw_white = fullfile(ruta_captura, ['WHITEREF_', nombre_muestra, '.raw']);
        if ~exist(hdr_white, 'file')
            hdr_white = fullfile(ruta_captura, ['WHITEREF_FX17_', nombre_muestra, '.hdr']);
            raw_white = fullfile(ruta_captura, ['WHITEREF_FX17_', nombre_muestra, '.raw']);
        end

        hdr_dark  = fullfile(ruta_captura, ['DARKREF_',  nombre_muestra, '.hdr']);
        raw_dark  = fullfile(ruta_captura, ['DARKREF_',  nombre_muestra, '.raw']);
        if ~exist(hdr_dark, 'file')
            hdr_dark = fullfile(ruta_captura, ['DARKREF_FX17_', nombre_muestra, '.hdr']);
            raw_dark = fullfile(ruta_captura, ['DARKREF_FX17_', nombre_muestra, '.raw']);
        end

        info  = read_info(hdr_muestra);
        cube  = read_cube(info,  raw_muestra);
        infoW = read_info(hdr_white);    white = read_cube(infoW, raw_white);
        infoD = read_info(hdr_dark);     dark  = read_cube(infoD, raw_dark);

        nrow = min([size(cube,1), size(white,1), size(dark,1)]);
        ncol = min([size(cube,2), size(white,2), size(dark,2)]);
        nb   = min([size(cube,3), size(white,3), size(dark,3)]);
        cube  = cube (1:nrow, 1:ncol, 1:nb);
        white = white(1:nrow, 1:ncol, 1:nb);
        dark  = dark (1:nrow, 1:ncol, 1:nb);

        cube_cal = (double(cube) - double(dark)) ./ ...
                   max(eps, (double(white) - double(dark)));
        cube_cal(~isfinite(cube_cal)) = 0;
        cube_cal = max(0, min(cube_cal, 1));
    catch ME
        fprintf('  ERROR loading data for sample %s. Skipping. (%s)\n', ...
                nombre_muestra, ME.message);
        close all;
        continue;
    end

    %% --- Pseudo-RGB from calibrated cube ---
    idxR = max(1, round(0.85*nb));
    idxG = max(1, round(0.55*nb));
    idxB = max(1, round(0.25*nb));
    rgb_pseudo = cat(3, cube_cal(:,:,idxR), cube_cal(:,:,idxG), cube_cal(:,:,idxB));
    rgb_pseudo = rgb_pseudo - min(rgb_pseudo(:));
    mx  = max(rgb_pseudo(:));
    if mx > 0, rgb_pseudo = rgb_pseudo ./ mx; end
    rgb_pseudo = imadjust(rgb_pseudo, stretchlim(rgb_pseudo, [0.01 0.99]));

    %% --- Base HSV segmentation on pseudo-RGB ---
    hsv        = rgb2hsv(rgb_pseudo);
    h_channel  = hsv(:,:,1);
    s_channel  = hsv(:,:,2);
    v_channel  = hsv(:,:,3);

    mask_hue   = (h_channel > HUE_RANGE(1)) & (h_channel < HUE_RANGE(2));
    thresh_s   = graythresh(s_channel);
    mask_sat   = (s_channel > thresh_s) & ...
                 (s_channel > MIN_SATURATION_THRESHOLD_HSV);
    thresh_v_o = graythresh(v_channel);
    mask_val_low  = (v_channel > thresh_v_o * 0.5);
    mask_val_high = (v_channel < MAX_VALUE_THRESHOLD_HSV);

    mask_uvas_base = mask_hue & mask_sat & mask_val_low & mask_val_high;

    %% --- Improved morphological cleanup (STEP 1) ---
    mask_uvas_base = bwareaopen(mask_uvas_base, MIN_BERRY_AREA);
    se_close = strel('disk', 3);
    mask_uvas_base = imclose(mask_uvas_base, se_close);
    mask_uvas_base = imfill(mask_uvas_base, 'holes');
    se_open = strel('disk', 1);
    mask_uvas_base = imopen(mask_uvas_base, se_open);
    se_dil = strel('disk', 1);
    mask_uvas_base = imdilate(mask_uvas_base, se_dil);

    %% --- FIGURE 1: base mask ---
    hfig1 = figure('Name', 'STEP 1: BASE MASK (HSV anti-rachis)', ...
                   'WindowState','maximized', ...
                   'Color','w');
    imshow(mask_uvas_base);
    title(sprintf('STEP 1: Base Mask (HSV Anti-Rachis) - Sample %s', ...
          nombre_muestra));
    drawnow;

    img_out_path1 = fullfile(carpeta_salida12, ...
        sprintf('PASO1_Mascara_AntiRaquis_%s.png', nombre_muestra));
    exportgraphics(hfig1, img_out_path1, 'Resolution', 300);

    fprintf('STEP 1 completed.\n');

    %% ===============================
    %  === STEP 2: SAM on pseudo-RGB ===
    %% ===============================
    fprintf('\n>> STEP 2: Running SAM (imsegsam) on pseudo-RGB...\n');

    I_sam = im2uint8(rgb_pseudo);
    masksSAM   = imsegsam(I_sam);
    Lsam_raw   = labelmatrix(masksSAM);
    numLabels_raw = max(Lsam_raw(:));
    fprintf('  SAM raw: %d distinct labels (unfiltered)\n', numLabels_raw);

    [rows, cols, ~] = size(I_sam);
    if ~isequal(size(mask_uvas_base), [rows cols])
        error('mask_uvas_base size does not match the SAM image size.');
    end

    %% --- FIGURE 2: raw SAM, white background ---
    hfig_raw = figure('Name','SAM segmentation (raw, unfiltered)', ...
                      'WindowState','maximized', ...
                      'Color','w');
    imshow(label2rgb(Lsam_raw,'jet','w','shuffle'));
    title(sprintf('SAM (raw) - Sample %s', nombre_muestra));
    drawnow;

    img_out_raw = fullfile(carpeta_salida12, ...
        sprintf('PASO2_SAM_raw_%s.png', nombre_muestra));
    exportgraphics(hfig_raw, img_out_raw, 'Resolution', 300);

    %% --- Restrict SAM to the HSV grape mask ---
    Lsam = Lsam_raw;
    Lsam(~mask_uvas_base) = 0;

    %% --- FIRST PASS: estimate typical berry area ---
    labels_inside = unique(Lsam(:));
    labels_inside(labels_inside == 0) = [];
    nlabels = numel(labels_inside);

    areas_raw = zeros(nlabels, 1);
    for i = 1:nlabels
        lab = labels_inside(i);
        areas_raw(i) = nnz(Lsam == lab);
    end

    if isempty(areas_raw)
        warning('No SAM regions were found inside the mask. Skipping sample.');
        close all;
        continue;
    end

    medianArea      = median(areas_raw);
    MIN_AREA_BERRY  = 0.25 * medianArea;
    MAX_AREA_SIMPLE = 1.8  * medianArea;

    fprintf('  Median SAM area inside mask: %.1f px\n', medianArea);
    fprintf('  Threshold MIN_AREA_BERRY  = %.1f px\n', MIN_AREA_BERRY);
    fprintf('  Threshold MAX_AREA_SIMPLE = %.1f px\n', MAX_AREA_SIMPLE);

    %% --- SECOND PASS: relabel + split large blobs (watershed) ---
    Lberries = zeros(rows, cols, 'uint16');
    areas    = [];
    nuevo_id = 0;

    for i = 1:nlabels
        lab      = labels_inside(i);
        mask_lab = (Lsam == lab);
        area_lab = nnz(mask_lab);
        if area_lab == 0
            continue;
        end

        if area_lab <= MAX_AREA_SIMPLE
            nuevo_id = nuevo_id + 1;
            Lberries(mask_lab) = nuevo_id;
            areas(nuevo_id,1)  = area_lab; %#ok<AGROW>
        else
            D = bwdist(~mask_lab);
            D = -D;
            D(~mask_lab) = Inf;
            Llocal = watershed(D);

            sub_labels = unique(Llocal(mask_lab));
            sub_labels(sub_labels == 0) = [];

            for s = 1:numel(sub_labels)
                sublab   = sub_labels(s);
                mask_sub = mask_lab & (Llocal == sublab);
                area_sub = nnz(mask_sub);
                if area_sub < MIN_AREA_BERRY
                    continue;
                end

                nuevo_id = nuevo_id + 1;
                Lberries(mask_sub) = nuevo_id;
                areas(nuevo_id,1)  = area_sub; %#ok<AGROW>
            end
        end
    end

    num_bayas = nuevo_id;
    if num_bayas == 0
        warning('No labeled berries were generated after blob splitting. Skipping sample.');
        close all;
        continue;
    end

    fprintf('  Final labeled berries (after blob splitting): %d\n', num_bayas);
    if num_bayas > 0
        fprintf('  Final areas: min = %.1f, max = %.1f, median = %.1f px\n', ...
            min(areas), max(areas), median(areas));
    end

    %% --- FIGURE 3: labeled berries, white background ---
    hfig2 = figure('Name','STEP 2: SAM regions inside the mask (refined)', ...
                   'WindowState','maximized', ...
                   'Color','w');
    if num_bayas > 0
        rgb_labels = label2rgb(Lberries, 'jet', 'w', 'shuffle');
    else
        rgb_labels = uint8(255*ones(rows, cols, 3));
    end
    imshow(rgb_labels);
    title(sprintf('STEP 2: %d labeled berries - Sample %s', ...
          num_bayas, nombre_muestra));
    drawnow;

    img_out_path2 = fullfile(carpeta_salida12, ...
        sprintf('PASO2_Bayas_SAM_%s.png', nombre_muestra));
    exportgraphics(hfig2, img_out_path2, 'Resolution', 300);

    %% --- FIGURE 4: pseudo-RGB only where berries exist + contours ---
    rgb_bayas = rgb_pseudo;
    mask_all  = (Lberries > 0);
    if any(mask_all(:))
        for c = 1:3
            chan = rgb_bayas(:,:,c);
            chan(~mask_all) = 1;
            rgb_bayas(:,:,c) = chan;
        end
    else
        rgb_bayas = ones(size(rgb_bayas));
    end

    hfig3 = figure('Name','STEP 2: SAM contours (refined) on pseudo-RGB', ...
                   'WindowState','maximized', ...
                   'Color','w');
    imshow(rgb_bayas);
    hold on;

    if num_bayas > 0
        for idb = 1:num_bayas
            mask_id = (Lberries == idb);
            if ~any(mask_id(:)), continue; end
            B = bwboundaries(mask_id, 'noholes');
            for kB = 1:numel(B)
                boundary = B{kB};
                plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 0.8);
            end
        end
    end

    title(sprintf('STEP 2: Contours of SAM regions (refined) - Sample %s', nombre_muestra));
    hold off; drawnow;

    img_out_path3 = fullfile(carpeta_salida12, ...
        sprintf('PASO2_Bayas_ContornosRGB_SAM_%s.png', nombre_muestra));
    exportgraphics(hfig3, img_out_path3, 'Resolution', 300);

    %% --- Save label matrix and areas (traceability) ---
    labels_file = fullfile(carpeta_salida12, ...
        sprintf('labels_SAM_%s.mat', nombre_muestra));
    L_SAM_berries = Lberries;
    areas_SAM     = areas;
    save(labels_file, 'L_SAM_berries', 'areas_SAM');

    fprintf('\nSTEP 1 (HSV) and STEP 2 (SAM) completed for %s.\n', nombre_muestra);
    fprintf('  - Base mask:            %s\n', img_out_path1);
    fprintf('  - SAM raw:              %s\n', img_out_raw);
    fprintf('  - Labeled berries:      %s\n', img_out_path2);
    fprintf('  - Contours on RGB:      %s\n', img_out_path3);
    fprintf('  - Labels saved:         %s\n', labels_file);

    %% ==========================
    %  === STEP 3: top 90 berries + ROI = FULL BERRY eroded by 2 px ===
    %% ==========================
    fprintf('\n>> STEP 3: ROI = full berry with %d px erosion (top 90 berries by area)...\n', ERODE_PIX);

    % --- Geometric measures for all berries ---
    stats = regionprops(L_SAM_berries,'Area');
    num_bayas = numel(stats);
    Area_px   = reshape([stats.Area], [], 1);

    fprintf('  (STEP 3) regions in L_SAM_berries: %d\n', num_bayas);

    % --- Select the 90 largest berries ---
    [~, idxSort] = sort(Area_px,'descend');
    Nkeep        = min(90, num_bayas);
    idxTop       = idxSort(1:Nkeep);
    fprintf('  %d largest berries will be used.\n', Nkeep);

    % --- Edge erosion ---
    se_erode = strel('disk', ERODE_PIX);

    % --- Containers for labels/spectra ---
    L_center       = zeros(rows, cols, 'uint16'); % labels (1..Nkeep)
    Spectra_center = zeros(Nkeep, nb);            % mean spectra

    % --- Efficient cube indexing (avoids 3D repmat) ---
    cube2 = reshape(cube_cal, [], nb);  % (rows*cols) x nb

    for ii = 1:Nkeep
        idBerry = idxTop(ii);

        berry_mask = (L_SAM_berries == idBerry);

        % ROI: berry mask eroded by 2 px
        roi_mask = imerode(berry_mask, se_erode);

        % If erosion removes the ROI, use the full berry
        if nnz(roi_mask) < 20
            roi_mask = berry_mask;
        end

        % Label ROI
        L_center(roi_mask) = ii;

        % Mean spectrum
        roi_idx = roi_mask(:);
        Spectra_center(ii,:) = mean(cube2(roi_idx, :), 1);
    end

    %% --- FX17 wavelengths: fixed list provided by the user ---
    wl = [ ...
        937.33, 944.25, 951.16, 958.08, 965.00, 971.92, 978.85, 985.77, ...
        992.70, 999.63, 1006.57, 1013.50, 1020.44, 1027.38, 1034.32, 1041.27, ...
        1048.21, 1055.16, 1062.12, 1069.07, 1076.03, 1082.98, 1089.94, 1096.91, ...
        1103.87, 1110.84, 1117.81, 1124.78, 1131.75, 1138.73, 1145.71, 1152.69, ...
        1159.67, 1166.66, 1173.64, 1180.63, 1187.63, 1194.62, 1201.62, 1208.62, ...
        1215.62, 1222.62, 1229.63, 1236.63, 1243.64, 1250.66, 1257.67, 1264.69, ...
        1271.71, 1278.73, 1285.75, 1292.78, 1299.80, 1306.83, 1313.87, 1320.90, ...
        1327.94, 1334.98, 1342.02, 1349.06, 1356.11, 1363.16, 1370.21, 1377.26, ...
        1384.31, 1391.37, 1398.43, 1405.49, 1412.56, 1419.62, 1426.69, 1433.76, ...
        1440.84, 1447.91, 1454.99, 1462.07, 1469.15, 1476.23, 1483.32, 1490.41, ...
        1497.50, 1504.59, 1511.69, 1518.79, 1525.89, 1532.99, 1540.09, 1547.20, ...
        1554.31, 1561.42, 1568.53, 1575.65, 1582.77, 1589.89, 1597.01, 1604.13, ...
        1611.26, 1618.39, 1625.52, 1632.66, 1639.79, 1646.93, 1654.07, 1661.21, ...
        1668.36, 1675.51, 1682.65, 1689.81, 1696.96, 1704.12, 1711.28, 1718.44 ...
    ];

    if numel(wl) ~= nb
        warning('Number of wavelengths (%d) differs from nb (%d). wl(1:nb) will be used.', numel(wl), nb);
        wl = wl(1:nb);
    end

    %% FIGURE 5 (90 ROIs): ROIs on pseudo-RGB without background
    rgb_fig5 = rgb_pseudo;
    for c = 1:3
        chan = rgb_fig5(:,:,c);
        chan(~mask_uvas_base) = 1;   % white background outside STEP 1 mask
        rgb_fig5(:,:,c) = chan;
    end

    hfig5 = figure('Name','STEP 3: ROIs (full berry eroded) on pseudo-RGB (background removed)', ...
                   'WindowState','maximized', 'Color','w');
    imshow(rgb_fig5);
    hold on;

    for ii = 1:Nkeep
        mask_id = (L_center == ii);
        if ~any(mask_id(:)), continue; end
        B = bwboundaries(mask_id, 'noholes');
        for kB = 1:numel(B)
            boundary = B{kB};
            plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 0.8);
        end
    end

    title(sprintf('STEP 3: %d ROIs = full berries (erosion %d px) - Sample %s', ...
        Nkeep, ERODE_PIX, nombre_muestra));
    hold off; drawnow;

    %% FIGURE 6 (90 ROIs): mean spectra (quick inspection only)
    hfig6 = figure('Name','STEP 3: Mean spectra of ROIs (full berry eroded)', ...
                   'WindowState','maximized', 'Color','w');
    plot(wl, Spectra_center.', 'LineWidth', 0.8);
    xlabel('Wavelength (nm)');
    ylabel('Mean reflectance');
    title(sprintf('Mean spectra - %d ROIs (full berry, eroded %d px)', Nkeep, ERODE_PIX));
    grid on; box on;
    drawnow;

    fprintf('STEP 3 completed (figures only, no table export).\n');

    %% ==========================
    %  === STEP 4: Export (90 vs 70 ROIs) ===
    %% ==========================
    fprintf('\n>> STEP 4: Exporting spectra and figures (90 vs 70)...\n');

    % Output folders:
    carpeta_90 = fullfile(ruta_principal, 'resultados_PASO_3_ROIs_centrales_90');
    carpeta_70 = fullfile(ruta_principal, 'resultados_PASO_3_ROIs_centrales_70');
    if ~exist(carpeta_90, 'dir'); mkdir(carpeta_90); end
    if ~exist(carpeta_70, 'dir'); mkdir(carpeta_70); end

    %% 4.1 – Spectral QC: distance to the median
    medSpec  = median(Spectra_center, 1);      % 1 × nb
    diffSpec = Spectra_center - medSpec;       % Nkeep × nb
    dist2    = sum(diffSpec.^2, 2);
    dist_med = sqrt(dist2);                    % Euclidean distance

    [~, idxByDist] = sort(dist_med, 'ascend');
    Nclean   = min(70, Nkeep);
    idxClean = idxByDist(1:Nclean);

    fprintf('  Spectral QC: keeping the %d ROIs closest to the median.\n', Nclean);

    %% 4.2 – Excel table of the 90 ROIs (table only) in carpeta_90
    SampleName = repmat(string(nombre_muestra), Nkeep, 1);
    ROI_ID     = (1:Nkeep).';

    wl_str  = cellstr(num2str(wl(:), '%0.2f'));
    wl_str  = strrep(wl_str, ' ', '');
    wl_str  = strrep(wl_str, '.', 'p');
    wl_lbls = strcat('nm_', wl_str);
    wl_lbls = matlab.lang.makeUniqueStrings(wl_lbls);

    Tspec90 = table(SampleName, ROI_ID, dist_med, ...
                    'VariableNames', {'SampleName','ROI_ID','Dist_toMedian'});
    for j = 1:nb
        Tspec90.(wl_lbls{j}) = Spectra_center(:,j);
    end

    file_xlsx90 = fullfile(carpeta_90, ...
        sprintf('Spectra90_CentralROIs_%s.xlsx', nombre_muestra));
    writetable(Tspec90, file_xlsx90, 'Sheet', 'Spectra90', 'WriteMode','overwritesheet');

    fprintf('  [90 ROIs] Spectra saved to: %s\n', file_xlsx90);

    %% 4.3 – Table + figures of the 70 "clean" ROIs in carpeta_70
    SampleName70 = SampleName(idxClean);
    ROI_ID70     = ROI_ID(idxClean);
    dist70       = dist_med(idxClean);
    Spectra70    = Spectra_center(idxClean,:);

    %% --- Sample metadata for Excel70 ---
    tokens = split(string(nombre_muestra), "_");

    % Camera
    if numel(tokens) >= 1 && strlength(tokens(1)) > 0
        Camara = tokens(1);
    else
        Camara = string(nombre_muestra(1:min(4,end)));
    end

    % Variety
    Variedad = "";
    if numel(tokens) >= 2
        Variedad = tokens(2);
        vUp = upper(strtrim(Variedad));
        if vUp == "GN"
            Variedad = "Garnacha tinta";
        elseif vUp == "GB"
            Variedad = "Garnacha blanca";
        end
    end

    % Date
    Fecha = "";
    if numel(tokens) >= 3
        tokFecha = tokens(3);
        if strlength(tokFecha) == 4 && all(isstrprop(char(tokFecha),'digit'))
            dd = str2double(extractBetween(tokFecha,1,2));
            mm = str2double(extractBetween(tokFecha,3,4));
            if isfinite(dd) && isfinite(mm)
                Fecha = string(sprintf('%02d/%02d', dd, mm));
            else
                Fecha = tokFecha;
            end
        else
            Fecha = tokFecha;
        end
    end

    % SampleIndex
    SampleIndex_val = NaN;
    if numel(tokens) >= 1
        SampleIndex_val = str2double(tokens(end));
    end

    % Design lookup
    Tratamiento_val = "";
    Bloque_val      = "";
    Riego_val       = "";
    Manejo_val      = "";

    if isfinite(SampleIndex_val)
        hit = design(design.SampleIndex == SampleIndex_val, :);
        if ~isempty(hit)
            Tratamiento_val = hit.Tratamiento(1);
            Bloque_val      = hit.Bloque(1);
            Riego_val       = hit.Riego(1);
            Manejo_val      = hit.Manejo_cultivo(1);
        end
    end

    Nclean_local = numel(idxClean);

    Camara70      = repmat(string(Camara),      Nclean_local, 1);
    Variedad70    = repmat(string(Variedad),    Nclean_local, 1);
    Fecha70       = repmat(string(Fecha),       Nclean_local, 1);
    SampleIndex70 = repmat(double(SampleIndex_val), Nclean_local, 1);

    Tratamiento70 = repmat(string(Tratamiento_val), Nclean_local, 1);
    Bloque70      = repmat(string(Bloque_val),      Nclean_local, 1);
    Riego70       = repmat(string(Riego_val),       Nclean_local, 1);
    Manejo70      = repmat(string(Manejo_val),      Nclean_local, 1);

    %% --- Excel table for 70 ROIs with metadata ---
    Tspec70 = table( ...
        Camara70, Variedad70, Fecha70, SampleIndex70, ...
        Tratamiento70, Bloque70, Riego70, Manejo70, ...
        SampleName70, ROI_ID70, dist70, ...
        'VariableNames', { ...
            'Camara','Variedad','Fecha','SampleIndex', ...
            'Tratamiento','Bloque','Riego','Manejo_cultivo', ...
            'SampleName','ROI_ID','Dist_toMedian'} );

    for j = 1:nb
        Tspec70.(wl_lbls{j}) = Spectra70(:,j);
    end

    file_xlsx70 = fullfile(carpeta_70, ...
        sprintf('Spectra70_CentralROIs_%s.xlsx', nombre_muestra));
    writetable(Tspec70, file_xlsx70, 'Sheet', 'Spectra70', 'WriteMode','overwritesheet');

    fprintf('  [70 ROIs] Spectra saved to: %s\n', file_xlsx70);

    % --- Figure of 70 clean ROIs only ---
    rgb_fig5_70 = rgb_fig5;
    hfig5_70 = figure('Name','STEP 3: 70 clean ROIs (image)', ...
                      'WindowState','maximized', 'Color','w');
    imshow(rgb_fig5_70);
    hold on;
    for kk = 1:Nclean
        idROI = idxClean(kk);
        mask_id = (L_center == idROI);
        if ~any(mask_id(:)), continue; end
        B = bwboundaries(mask_id, 'noholes');
        for kB = 1:numel(B)
            boundary = B{kB};
            plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 0.8);
        end
    end
    title(sprintf('STEP 3: %d clean ROIs (full berry, eroded %d px) - Sample %s', ...
        Nclean, ERODE_PIX, nombre_muestra));
    hold off; drawnow;

    file_fig5_70 = fullfile(carpeta_70, ...
        sprintf('PASO3_ROIsCentrales70_%s.fig', nombre_muestra));
    file_jpg5_70 = fullfile(carpeta_70, ...
        sprintf('PASO3_ROIsCentrales70_%s.jpg', nombre_muestra));
    savefig(hfig5_70, file_fig5_70);
    exportgraphics(hfig5_70, file_jpg5_70, 'Resolution', 300);

    % --- Figure of 70 clean spectra only ---
    hfig6_70 = figure('Name','STEP 3: Mean spectra of 70 clean ROIs', ...
                      'WindowState','maximized', 'Color','w');
    plot(wl, Spectra70.', 'LineWidth', 0.8);
    xlabel('Wavelength (nm)');
    ylabel('Mean reflectance');
    title(sprintf('Mean spectra - %d clean ROIs (full berry, eroded %d px)', Nclean, ERODE_PIX));
    grid on; box on; drawnow;

    file_fig6_70 = fullfile(carpeta_70, ...
        sprintf('PASO3_MeanSpectra70_%s.fig', nombre_muestra));
    file_jpg6_70 = fullfile(carpeta_70, ...
        sprintf('PASO3_MeanSpectra70_%s.jpg', nombre_muestra));
    savefig(hfig6_70, file_fig6_70);
    exportgraphics(hfig6_70, file_jpg6_70, 'Resolution', 300);

    fprintf('  [70 ROIs] Figures saved in:\n    %s\n    %s\n    %s\n    %s\n', ...
        file_fig5_70, file_jpg5_70, file_fig6_70, file_jpg6_70);

    %% 4.4 – Close all figures for this sample
    close all;

    fprintf('STEP 4 completed for %s.\n', nombre_muestra);
    fprintf('>> End of sample.\n');

end % end sample loop

fprintf('\n>> End of script: all FX17_ samples processed.\n');

%% ==========================================================
%  === STEP 5: COMPILE 70-ROI Excel files + sample means
%% ==========================================================
fprintf('\n>> STEP 5: Compiling 70-ROI Excel files and calculating sample means...\n');

carpeta_70 = fullfile(ruta_principal, 'resultados_PASO_3_ROIs_centrales_70');
files70 = dir(fullfile(carpeta_70, 'Spectra70_CentralROIs_*.xlsx'));

if isempty(files70)
    warning('No Spectra70_CentralROIs_*.xlsx files were found in %s. Results will not be compiled.', carpeta_70);
    return;
end

Tall = table();
Tmeans = table();

for iF = 1:numel(files70)
    fpath = fullfile(files70(iF).folder, files70(iF).name);
    fprintf('  Reading: %s\n', files70(iF).name);

    try
        T = readtable(fpath, 'Sheet', 'Spectra70', 'TextType', 'string');
    catch
        T = readtable(fpath, 'Sheet', 'Spectra70');
    end

    % Accumulate all rows
    if isempty(Tall)
        Tall = T;
    else
        vAll = Tall.Properties.VariableNames;
        vT   = T.Properties.VariableNames;

        if ~isequal(vAll, vT)
            missingInT = setdiff(vAll, vT);
            for m = 1:numel(missingInT)
                TallVar = missingInT{m};
                T.(TallVar) = missing(height(T), 1);
            end

            missingInAll = setdiff(vT, vAll);
            for m = 1:numel(missingInAll)
                Tvar = missingInAll{m};
                Tall.(Tvar) = missing(height(Tall), 1);
            end

            T = T(:, Tall.Properties.VariableNames);
        end

        Tall = [Tall; T]; %#ok<AGROW>
    end

    % --- Calculate sample mean ---
    varNames = T.Properties.VariableNames;
    isSpec = startsWith(varNames, 'nm_');

    metaVarsWanted = {'Camara','Variedad','Fecha','SampleIndex','Tratamiento','Bloque','Riego','Manejo_cultivo','SampleName'};
    metaVars = intersect(metaVarsWanted, varNames, 'stable');

    Tm = table();

    for mv = 1:numel(metaVars)
        v = metaVars{mv};
        Tm.(v) = T.(v)(1);
    end

    Tm.N_ROIs = height(T);

    if any(strcmp(varNames, 'Dist_toMedian'))
        try
            Tm.Mean_Dist_toMedian = mean(T.Dist_toMedian, 'omitnan');
        catch
            Tm.Mean_Dist_toMedian = mean(double(T.Dist_toMedian), 'omitnan');
        end
    end

    if any(isSpec)
        X = T{:, isSpec};
        X = double(X);
        mu = mean(X, 1, 'omitnan');
        specNames = varNames(isSpec);
        for j = 1:numel(specNames)
            Tm.(specNames{j}) = mu(j);
        end
    end

    if isempty(Tmeans)
        Tmeans = Tm;
    else
        Tmeans = [Tmeans; Tm]; %#ok<AGROW>
    end
end

% Save compiled table
out_compiled = fullfile(carpeta_70, 'Spectra70_ALL_Samples_compiled.xlsx');
writetable(Tall, out_compiled, 'Sheet', 'All_70ROIs', 'WriteMode', 'overwritesheet');
fprintf('  [OK] Compiled table (all rows) saved to:\n    %s\n', out_compiled);

% Save sample means
out_means = fullfile(carpeta_70, 'Spectra70_MEANS_bySample.xlsx');
writetable(Tmeans, out_means, 'Sheet', 'Means', 'WriteMode', 'overwritesheet');
fprintf('  [OK] Sample means saved to:\n    %s\n', out_means);

fprintf('>> STEP 5 completed.\n');

end