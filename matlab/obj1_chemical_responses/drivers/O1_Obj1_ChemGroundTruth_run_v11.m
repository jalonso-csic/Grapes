function O1_Obj1_ChemGroundTruth_run_v11(inputXlsx, sheetName)
% O1_Obj1_ChemGroundTruth_run_v11 (audit-ready)
% ------------------------------------------------------------
% Objective 1: Chemical ground-truth statistics for the grape dataset
%
% This script:
%   (i) reads the master chemistry table (one row per experimental unit),
%   (ii) harmonises design factors (Date, Irrigation, Management, Block),
%   (iii) fits linear mixed-effects models (LMMs) per endpoint,
%   (iv) computes planned Irrigation–Rainfed contrasts by date,
%   (v) exports manuscript-oriented tables plus publication-quality figures.
%
% OUTPUTS (written under a folder named "Objetivo_1"):
%   - O1_Stats_Summary.xlsx     (tidy ANOVA + contrasts + means)
%   - O1_Table1_Manuscript.xlsx (Table 1-ready p-values)
%   - Fig*.png (300 dpi) + Fig*.fig
%   - O1_log.txt
%
% USAGE:
%   O1_Obj1_ChemGroundTruth_run_v11('Matriz_CHEM_HSI_MASTER_96.xlsx','Matriz');

    % ---------------------------
    % 0) Robust inputs
    % ---------------------------
    if nargin < 1 || isempty(inputXlsx)
        inputXlsx = 'Matriz_CHEM_HSI_MASTER_96.xlsx';
    end
    if exist(inputXlsx,'file') ~= 2
        error('Input file not found: %s', inputXlsx);
    end

    if nargin < 2 || isempty(sheetName)
        try
            sh = sheetnames(inputXlsx);
            sheetName = sh{1};
        catch
            sheetName = 'Matriz';
        end
    end

    ts = datestr(now,'yyyymmdd_HHMMSS');

    % Parent output folder requested by the user:
    baseDir = fileparts(inputXlsx);
    if isempty(baseDir); baseDir = pwd; end
    objDir  = fullfile(baseDir, 'Objetivo_1');
    if ~exist(objDir,'dir'); mkdir(objDir); end

    outDir = fullfile(objDir, ['RESULTS_Obj1_' ts]);
    if ~exist(outDir,'dir'); mkdir(outDir); end

    % ---------------------------
    % 1) Log
    % ---------------------------
    logFile = fullfile(outDir,'O1_log.txt');
    diary(logFile); diary on;
    fprintf('=== Objective 1 | Chemical Ground Truth (v11) ===\n');
    fprintf('Input: %s\n', inputXlsx);
    fprintf('Sheet: %s\n', sheetName);
    fprintf('Timestamp: %s\n\n', ts);

    % ---------------------------
    % 2) Read table (preserve headers)
    % ---------------------------
    Traw = readtable_preserve(inputXlsx, sheetName);

    % ---------------------------
    % 3) Resolve design factors
    % ---------------------------
    V.Date  = resolveVar(Traw, {'SamplingDate','Fecha','Fecha_ddmm','Sampling_Event_Date'});
    V.Irrig = resolveVar(Traw, {'IrrigationRegime','Riego','Irrigation'});
    V.Mgmt  = resolveVar(Traw, {'VineyardFloorManagement','Manejo_cultivo','Management'});
    V.Trt   = resolveVar(Traw, {'Treatment','Tratamiento'});
    V.Block = resolveVar(Traw, {'Block','BloqueRoman','Bloque','Replicate'});

    mustHave = {'Date','Irrig','Mgmt','Trt','Block'};
    for k=1:numel(mustHave)
        if strcmp(V.(mustHave{k}), "MISSING")
            error('Missing required column: %s (check Excel headers).', mustHave{k});
        end
    end

    % Standardise working factors
    T = table();
    T.DateRaw  = Traw.(V.Date);
    T.IrrigRaw = Traw.(V.Irrig);
    T.MgmtRaw  = Traw.(V.Mgmt);
    T.TrtRaw   = Traw.(V.Trt);
    T.BlockRaw = Traw.(V.Block);

    T.Date  = makeDateCategorical(T.DateRaw, 2025);     % dd/MM ordered categorical
    T.Irrig = makeIrrigCategorical(T.IrrigRaw);         % Irrigation / Rainfed
    T.Mgmt  = makeMgmtCategorical(T.MgmtRaw);           % standardised levels where possible
    T.Trt   = categorical(string(T.TrtRaw));
    T.Block = makeBlockCategorical(T.BlockRaw);         % I–IV (when numeric) or preserved strings

    fprintf('\n[DIAG] Date categories:\n');
    disp(categories(T.Date));
    disp(countcats(T.Date));
    fprintf('[DIAG] Irrigation categories:\n');
    disp(categories(T.Irrig));
    disp(countcats(T.Irrig));

    % ---------------------------
    % 4) Endpoints
    % ---------------------------
    EP = endpointList();
    EP = resolveAllEndpoints(Traw, EP);

    fprintf('\n[DIAG] Resolved endpoints:\n');
    for i=1:numel(EP)
        fprintf('  %-6s -> %s\n', EP(i).code, EP(i).var);
    end
    fprintf('\n');

    % ---------------------------
    % 5) Output workbooks
    % ---------------------------
    statsXlsx  = fullfile(outDir,'O1_Stats_Summary.xlsx');
    table1Xlsx = fullfile(outDir,'O1_Table1_Manuscript.xlsx');

    ANOVA_Global_Irrig = table();
    ANOVA_Irrig_Mgmt   = table();
    ANOVA_Rain_Mgmt    = table();
    Means_Irrig        = table();
    Means_Mgmt_Irrig   = table();
    Means_Mgmt_Rain    = table();
    Contr_Irrig_byDate = table();

    warning('off','stats:LinearMixedModel:RankDeficientDesign');
    warning('off','stats:LinearMixedModel:IgnoreCovariancePattern');

    % ---------------------------
    % 6) Loop over endpoints
    % ---------------------------
    for i = 1:numel(EP)
        code = EP(i).code;
        yvar = EP(i).var;
        unit = EP(i).unit;
        [plotLabel, plotUnit] = getEnglishLabel(code);

        fprintf('Processing: %-6s -> %s [%s]\n', code, plotLabel, yvar);

        if ~ismember(yvar, string(Traw.Properties.VariableNames))
            fprintf('  [skip] Missing variable: %s\n', yvar);
            continue;
        end

        D = T;
        D.y = Traw.(yvar);

        % Remove missing y and also missing predictors
        D = rmmissing(D,'DataVariables','y');
        D = rmmissing(D);

        % Drop unused categories after filtering (important for stable design matrices)
        D.Date  = removecats(D.Date);
        D.Irrig = removecats(D.Irrig);
        D.Mgmt  = removecats(D.Mgmt);
        D.Block = removecats(D.Block);

        if height(D) < 12
            fprintf('  [skip] Too few rows after filtering.\n');
            continue;
        end

        if numel(categories(D.Date)) < 2 || numel(categories(D.Irrig)) < 2 || numel(categories(D.Block)) < 2
            fprintf('  [skip] Not enough factor levels (Date/Irrigation/Block).\n');
            continue;
        end

        % -------- Means: Date × Irrigation
        mI = meansSE(D, {'Date','Irrig'}, 'y');
        mI.Endpoint = repmat(string(code), height(mI), 1);
        Means_Irrig = [Means_Irrig; mI]; %#ok<AGROW>

        % -------- Subsets
        Di = D(D.Irrig=="Irrigation",:);
        Dr = D(D.Irrig=="Rainfed",:);

        % Drop unused categories after subsetting
        if ~isempty(Di)
            Di.Date  = removecats(Di.Date);
            Di.Mgmt  = removecats(Di.Mgmt);
            Di.Block = removecats(Di.Block);
        end
        if ~isempty(Dr)
            Dr.Date  = removecats(Dr.Date);
            Dr.Mgmt  = removecats(Dr.Mgmt);
            Dr.Block = removecats(Dr.Block);
        end

        % -------- Means: Date × Management within irrigated
        if height(Di) > 0 && numel(categories(Di.Mgmt)) >= 2
            mMi = meansSE(Di, {'Date','Mgmt'}, 'y');
            mMi.Endpoint = repmat(string(code), height(mMi), 1);
            Means_Mgmt_Irrig = [Means_Mgmt_Irrig; mMi]; %#ok<AGROW>
        else
            mMi = table();
        end

        % -------- Means: Date × Management within rainfed
        if height(Dr) > 0 && numel(categories(Dr.Mgmt)) >= 2
            mMr = meansSE(Dr, {'Date','Mgmt'}, 'y');
            mMr.Endpoint = repmat(string(code), height(mMr), 1);
            Means_Mgmt_Rain = [Means_Mgmt_Rain; mMr]; %#ok<AGROW>
        else
            mMr = table();
        end

        % -------- LMM: Global (Date×Irrigation + (1|Block))
        try
            lme_g = fitlme(D, 'y ~ Date*Irrig + (1|Block)', 'FitMethod','REML');
            a_g = anovaToTidy(toTable(anova(lme_g,'DFMethod','satterthwaite')), code, unit, "Global_DatexIrrig");
            a_g.Term = harmoniseTerms(a_g.Term);
            ANOVA_Global_Irrig = [ANOVA_Global_Irrig; a_g]; %#ok<AGROW>
        catch ME
            fprintf('  [warn] Global LMM failed: %s\n', ME.message);
        end

        % -------- LMM: Management within irrigated
        if height(Di) > 12 && numel(categories(Di.Mgmt)) >= 2 && numel(categories(Di.Date)) >= 2
            [okI, a_i] = fitlme_anova_try(Di, code, unit, "Irrigated", true);
            if okI
                a_i.Term = harmoniseTerms(a_i.Term);
                ANOVA_Irrig_Mgmt = [ANOVA_Irrig_Mgmt; a_i]; %#ok<AGROW>
            end
        end

        % -------- LMM: Management within rainfed
        if height(Dr) > 12 && numel(categories(Dr.Mgmt)) >= 2 && numel(categories(Dr.Date)) >= 2
            [okR, a_r] = fitlme_anova_try(Dr, code, unit, "Rainfed", true);
            if okR
                a_r.Term = harmoniseTerms(a_r.Term);
                ANOVA_Rain_Mgmt = [ANOVA_Rain_Mgmt; a_r]; %#ok<AGROW>
            end
        end

        % -------- Contrasts: Irrigation vs Rainfed by date
        try
            C1 = contrast_2level_byDate(D, 'Date', 'y', 'Irrig', 'Block', "Irrigation", "Rainfed");
            C1.Endpoint = repmat(string(code), height(C1), 1);
            Contr_Irrig_byDate = [Contr_Irrig_byDate; C1]; %#ok<AGROW>
        catch ME
            fprintf('  [warn] Contrast-by-date failed: %s\n', ME.message);
        end

        % ---------------------------
        % Figures (publication-oriented)
        % ---------------------------
        fig1 = plotMeansSE_Pub(mI, 'Date', 'Irrig', plotLabel, plotUnit);
        exportgraphics(fig1, fullfile(outDir, sprintf('Fig1_Irrigation_%s.png', code)), ...
            'Resolution',300,'BackgroundColor','w');
        savefig(fig1, fullfile(outDir, sprintf('Fig1_Irrigation_%s.fig', code)));
        close(fig1);

        if ~isempty(mMi) && height(mMi) > 0
            fig2 = plotMeansSE_Pub(mMi, 'Date', 'Mgmt', plotLabel, plotUnit);
            exportgraphics(fig2, fullfile(outDir, sprintf('Fig2_Management_Irrigated_%s.png', code)), ...
                'Resolution',300,'BackgroundColor','w');
            savefig(fig2, fullfile(outDir, sprintf('Fig2_Management_Irrigated_%s.fig', code)));
            close(fig2);
        end

        if ~isempty(mMr) && height(mMr) > 0
            fig3 = plotMeansSE_Pub(mMr, 'Date', 'Mgmt', plotLabel, plotUnit);
            exportgraphics(fig3, fullfile(outDir, sprintf('Fig3_Management_Rainfed_%s.png', code)), ...
                'Resolution',300,'BackgroundColor','w');
            savefig(fig3, fullfile(outDir, sprintf('Fig3_Management_Rainfed_%s.fig', code)));
            close(fig3);
        end
    end

    warning('on','stats:LinearMixedModel:RankDeficientDesign');
    warning('on','stats:LinearMixedModel:IgnoreCovariancePattern');

    % ---------------------------
    % 7) Write stats workbook
    % ---------------------------
    fprintf('\nWriting Excel summary...\n');
    writetable_safe(ANOVA_Global_Irrig, statsXlsx, 'ANOVA_Global_DatexIrrig');
    writetable_safe(ANOVA_Irrig_Mgmt,   statsXlsx, 'ANOVA_Irrigated_Mgmt');
    writetable_safe(ANOVA_Rain_Mgmt,    statsXlsx, 'ANOVA_Rainfed_Mgmt');
    writetable_safe(Contr_Irrig_byDate, statsXlsx, 'Contr_Irrig_vs_Rainfed');
    writetable_safe(Means_Irrig,        statsXlsx, 'Means_DatexIrrig');
    writetable_safe(Means_Mgmt_Irrig,   statsXlsx, 'Means_DatexMgmt_Irrig');
    writetable_safe(Means_Mgmt_Rain,    statsXlsx, 'Means_DatexMgmt_Rain');

    % ---------------------------
    % 8) Manuscript-ready Table 1
    % ---------------------------
    T1 = buildTable1Manuscript(EP, ANOVA_Global_Irrig, ANOVA_Irrig_Mgmt, ANOVA_Rain_Mgmt);
    try
        writetable(T1, table1Xlsx, 'Sheet','Table1', 'WriteMode','overwritesheet');
    catch
        writetable(T1, table1Xlsx, 'Sheet','Table1');
    end

    fprintf('\n[OK] Done. Results in: %s\n', outDir);
    diary off;
end

% =========================================================
%  ENDPOINTS
% =========================================================
function EP = endpointList()
    EP = repmat(struct('code',"",'label',"",'unit',"",'domain',"",'candidates',{{}},'var',""), 0, 1);
    EP(end+1) = ep('BW',    'Berry weight (50 berries)',           'g',      'Tech', {'BerryWeight50_g','Peso_50_bayas_g'});
    EP(end+1) = ep('Brix',  'Total soluble solids',                '°Brix',  'Tech', {'TSS_Brix','Sólidos solubles totales (ºBrix)','Brix'});
    EP(end+1) = ep('pH',    'pH',                                  '',       'Tech', {'pH'});
    EP(end+1) = ep('AT',    'Titratable acidity',                  'g/L',    'Tech', {'TA_gL','Acidez titulable (g/L)','AT'});
    EP(end+1) = ep('Malic', 'Malic acid',                          'g/L',    'Tech', {'MalicAcid_gL','A.Málico (g/L)','A.Malico (g/L)'});
    EP(end+1) = ep('Tart',  'Tartaric acid',                       'g/L',    'Tech', {'TartaricAcid_gL','A.Tartárico (g/L)','A.Tartarico (g/L)'});
    EP(end+1) = ep('K',     'Potassium',                           'g/L',    'Tech', {'Potassium_gL','POTASIO (g/L)','Potasio'});
    EP(end+1) = ep('NOPA',  'YAN (alpha-amino N)',                 'mg N/L', 'Tech', {'AlphaAminoN_mgL','Nitrógeno α-amínico (NOPA, mg/L)','NOPA'});
    EP(end+1) = ep('NH4',   'Ammonium',                            'mg N/L', 'Tech', {'Ammonium_mgL','Amonio (NH4+, mg/L)','Amonio'});
    EP(end+1) = ep('PPT',   'Total phenolic potential',            'mg/kg',  'Phen', {'TotalPhenolicPotential_mgkg','PPT_mgkg_mean','PPT'});
    EP(end+1) = ep('PPEx',  'Extractable phenolic potential',      'mg/kg',  'Phen', {'ExtractablePhenolicPotential_mgkg','PPEx_mgkg_mean','PPEx'});
    EP(end+1) = ep('AnT',   'Total anthocyanins',                  'mg/kg',  'Phen', {'TotalAnthocyanins_mgkg','AnT_mgkg_mean','AnT'});
    EP(end+1) = ep('AnEx',  'Extractable anthocyanins',            'mg/kg',  'Phen', {'ExtractableAnthocyanins_mgkg','AnEx_mgkg_mean','AnEx'});
end

function s = ep(code,label,unit,domain,candidates)
    s = struct('code',string(code),'label',string(label),'unit',string(unit), ...
        'domain',string(domain),'candidates',{candidates},'var',string(""));
end

function [label, unit] = getEnglishLabel(code)
    switch char(code)
        case 'BW',   label='Berry weight (50 berries)'; unit='g';
        case 'Brix', label='Total soluble solids'; unit='°Brix';
        case 'pH',   label='pH'; unit='';
        case 'AT',   label='Titratable acidity'; unit='g/L tartaric eq';
        case 'Malic',label='Malic acid'; unit='g/L';
        case 'Tart', label='Tartaric acid'; unit='g/L';
        case 'K',    label='Potassium'; unit='g/L';
        case 'NOPA', label='YAN (alpha-amino N)'; unit='mg N/L';
        case 'NH4',  label='Ammonium'; unit='mg N/L';
        case 'PPT',  label='Total phenolic potential'; unit='mg/kg';
        case 'PPEx', label='Extractable phenolic potential'; unit='mg/kg';
        case 'AnT',  label='Total anthocyanins'; unit='mg/kg';
        case 'AnEx', label='Extractable anthocyanins'; unit='mg/kg';
        otherwise,   label=char(code); unit='';
    end
end

% =========================================================
%  LMM helper with optional fallback
% =========================================================
function [ok, a] = fitlme_anova_try(Dsub, endpointCode, unit, tagPrefix, allowFallback)
    ok = false; a = table();

    % Primary: interaction
    try
        lme = fitlme(Dsub, 'y ~ Date*Mgmt + (1|Block)', 'FitMethod','REML');
        a = anovaToTidy(toTable(anova(lme,'DFMethod','satterthwaite')), endpointCode, unit, tagPrefix + "_DatexMgmt");
        ok = true;
        return;
    catch ME
        fprintf('  [warn] %s management LMM (Date*Management) failed: %s\n', tagPrefix, ME.message);
    end

    if ~allowFallback, return; end

    % Fallback: main effects only
    try
        lme = fitlme(Dsub, 'y ~ Date + Mgmt + (1|Block)', 'FitMethod','REML');
        a = anovaToTidy(toTable(anova(lme,'DFMethod','satterthwaite')), endpointCode, unit, tagPrefix + "_DatePlusMgmt");
        ok = true;
        fprintf('  [info] %s management LMM used fallback: Date + Management\n', tagPrefix);
    catch ME
        fprintf('  [warn] %s management LMM fallback (Date+Management) failed: %s\n', tagPrefix, ME.message);
    end
end

% =========================================================
%  Variable resolution / reading
% =========================================================
function v = resolveVar(T, candidates)
    vnames = string(T.Properties.VariableNames);
    for k=1:numel(candidates)
        if any(vnames == string(candidates{k}))
            v = string(candidates{k});
            return;
        end
    end
    vnorm = normalizeTokens(vnames);
    cnorm = normalizeTokens(string(candidates));
    for k=1:numel(cnorm)
        idx = find(contains(vnorm, cnorm(k)), 1);
        if ~isempty(idx)
            v = vnames(idx);
            return;
        end
    end
    v = "MISSING";
end

function T = readtable_preserve(xlsx, sheetName)
    opts = detectImportOptions(xlsx, 'Sheet', sheetName);
    opts.VariableNamingRule = 'preserve';
    T = readtable(xlsx, opts);
end

function x = normalizeTokens(x)
    x = lower(string(x));
    x = replace(x, "α","a");
    x = replace(x, ["á","à","ä","â"], "a");
    x = replace(x, ["é","è","ë","ê"], "e");
    x = replace(x, ["í","ì","ï","î"], "i");
    x = replace(x, ["ó","ò","ö","ô"], "o");
    x = replace(x, ["ú","ù","ü","û"], "u");
    x = replace(x, "ñ","n");
    x = regexprep(x, '[^a-z0-9]+','');
end

% =========================================================
%  Factor builders
% =========================================================
function C = makeDateCategorical(x, year0)
    if isdatetime(x)
        dt = x;
    elseif isnumeric(x)
        if all(x(~isnan(x)) > 20000)
            dt = datetime(x, 'ConvertFrom','excel');
        else
            error('Numeric SamplingDate not recognised (not Excel serial).');
        end
    else
        s = string(x); s = strtrim(s);

        dt = datetime(s, 'InputFormat','yyyy-MM-dd');
        if any(isnat(dt)), dt = datetime(s, 'InputFormat','dd/MM/yyyy'); end
        if any(isnat(dt)), dt = datetime(s, 'InputFormat','dd-MM-yyyy'); end

        if any(isnat(dt))
            dt = NaT(size(s));
            for i=1:numel(s)
                if strlength(s(i))==0, continue; end
                parts = split(s(i), ["/","-"]);
                if numel(parts)==2
                    dd = str2double(parts(1)); mm = str2double(parts(2));
                    if ~isnan(dd) && ~isnan(mm)
                        dt(i) = datetime(year0, mm, dd);
                    end
                elseif numel(parts)==3
                    dd = str2double(parts(1)); mm = str2double(parts(2)); yy = str2double(parts(3));
                    if ~isnan(dd) && ~isnan(mm) && ~isnan(yy)
                        dt(i) = datetime(yy, mm, dd);
                    end
                end
            end
        end
    end

    if any(isnat(dt))
        error('SamplingDate contains unparsable values -> would create <undefined>.');
    end

    u = unique(dt); u = sort(u);
    cats = compose('%02d/%02d', day(u), month(u));
    vals = compose('%02d/%02d', day(dt), month(dt));
    C = categorical(vals, cats, cats, 'Ordinal', true);
end

function C = makeIrrigCategorical(x)
    s = lower(string(x));
    s = strtrim(s);
    out = strings(size(s));
    out(contains(s,"riego") | contains(s,"irrig")) = "Irrigation";
    out(contains(s,"secano") | contains(s,"rain")) = "Rainfed";
    out(out=="") = s(out=="");
    C = categorical(out);
    if all(ismember(["Rainfed","Irrigation"], categories(C)))
        C = reordercats(C, ["Rainfed","Irrigation"]);
    end
    C = removecats(C);
end

function C = makeMgmtCategorical(x)
    s = lower(string(x));
    s = strtrim(s);
    s(contains(s,'open'))       = "open vegetation";
    s(contains(s,'tillage'))    = "tillage";
    s(contains(s,'alternate'))  = "alternate cover crop";
    s(contains(s,'double'))     = "double cover crop";
    s(contains(s,'vegetacion')) = "open vegetation";
    s(contains(s,'laboreo'))    = "tillage";
    s(contains(s,'alterna'))    = "alternate cover crop";
    s(contains(s,'doble'))      = "double cover crop";
    C = categorical(s);
    pref = ["open vegetation","tillage","alternate cover crop","double cover crop"];
    present = intersect(pref, string(categories(C)), 'stable');
    others  = setdiff(string(categories(C)), present, 'stable');
    if ~isempty(present)
        C = reordercats(C, [present(:); others(:)]);
    end
    C = removecats(C);
end

function C = makeBlockCategorical(x)
    if isnumeric(x)
        xi = x(:);
        rom = strings(size(xi));
        rom(xi==1)="I"; rom(xi==2)="II"; rom(xi==3)="III"; rom(xi==4)="IV";
        rom(rom=="") = string(xi(rom==""));
        C = categorical(rom);
    else
        C = categorical(string(x));
    end
    C = removecats(C);
end

% =========================================================
%  Endpoint resolution
% =========================================================
function EP = resolveAllEndpoints(T, EP)
    vnames = string(T.Properties.VariableNames);
    vnorm  = normalizeTokens(vnames);
    for i = 1:numel(EP)
        cands = string(EP(i).candidates);
        hit = "";
        for k = 1:numel(cands)
            if any(vnames == cands(k)); hit = cands(k); break; end
        end
        if hit == ""
            cnorm = normalizeTokens(cands);
            for k = 1:numel(cnorm)
                idx = find(contains(vnorm, cnorm(k)), 1);
                if ~isempty(idx); hit = vnames(idx); break; end
            end
        end
        if hit ~= ""; EP(i).var = hit; else; EP(i).var = "MISSING"; end
    end
    found = [EP.var]; EP(found=="MISSING") = [];
end

% =========================================================
%  Stats helpers
% =========================================================
function G = meansSE(T, groupVars, yName)
    y = T.(yName);
    [gid, Glevels] = findgroups(T(:, groupVars));
    mu = splitapply(@(x) mean(x,'omitnan'), y, gid);
    sd = splitapply(@(x) std(x,0,'omitnan'), y, gid);
    n  = splitapply(@(x) sum(~ismissing(x)), y, gid);
    se = sd ./ sqrt(max(1,n));
    G = Glevels;
    G.mean = mu; G.sd = sd; G.n = n; G.se = se;
end

function tidy = anovaToTidy(a, endpointCode, unit, modelTag)
    v = string(a.Properties.VariableNames);
    termCol = find(v=="Term" | v=="Terms" | v=="Source", 1);
    pCol    = find(v=="pValue" | v=="Prob_F" | v=="ProbF" | v=="p", 1);
    FCol    = find(v=="FStat" | v=="F", 1);
    if isempty(termCol), termCol = 1; end
    if isempty(pCol), p = nan(height(a),1); else, p = a{:,pCol}; end
    if isempty(FCol), F = nan(height(a),1); else, F = a{:,FCol}; end
    tidy = table();
    tidy.Endpoint = repmat(string(endpointCode), height(a), 1);
    tidy.Unit     = repmat(string(unit), height(a), 1);
    tidy.Model    = repmat(string(modelTag), height(a), 1);
    tidy.Term     = string(a{:,termCol});
    tidy.F        = F;
    tidy.pValue   = p;
end

function A = toTable(x)
    if istable(x), A = x; return; end
    if isa(x,'dataset'), A = dataset2table(x); return; end
    A = struct2table(struct(x));
end

function t = harmoniseTerms(t)
    s = string(t);
    s = replace(s, "Irrig", "Irrigation");
    s = replace(s, "Mgmt",  "Management");
    s = replace(s, "Date:Irrigation", "Date×Irrigation");
    s = replace(s, "Date:Irrig",      "Date×Irrigation");
    s = replace(s, "Date:Management", "Date×Management");
    s = replace(s, "Date:Mgmt",       "Date×Management");
    t = s;
end

function C = contrast_2level_byDate(T, dateVar, yName, groupVar, blockVar, levelA, levelB)
    dates = categories(T.(dateVar));
    C = table();
    for i = 1:numel(dates)
        d = dates{i};
        Td = T(T.(dateVar)==d, :);

        g = Td.(groupVar);
        present = unique(string(g(~ismissing(g))));
        if ~all(ismember([string(levelA), string(levelB)], present)), continue; end

        g = removecats(g);
        other = setdiff(categories(g), {char(levelB)}, 'stable');
        g = reordercats(g, [{char(levelB)}; other(:)]);
        Td.(groupVar) = g;

        f = sprintf('%s ~ %s + (1|%s)', yName, groupVar, blockVar);
        lme = fitlme(Td, f, 'FitMethod','REML');

        coefTab = toTable(lme.Coefficients);
        cNames = string(coefTab.Name);
        row = find(contains(cNames, string(levelA)) & contains(cNames, string(groupVar)), 1);

        if ~isempty(row)
            est = coefTab.Estimate(row); se=coefTab.SE(row); df=coefTab.DF(row);
            tstat=coefTab.tStat(row); p=coefTab.pValue(row);
        else
            est=nan; se=nan; df=nan; tstat=nan; p=nan;
        end

        tcrit = tinv(0.975, max(1,df));
        rowT = table(string(d), string(levelA), string(levelB), est, se, df, tstat, p, ...
            est-tcrit*se, est+tcrit*se, ...
            'VariableNames', {'Date','LevelA','LevelB','Estimate','SE','DF','tStat','pValue','CI95_Lo','CI95_Hi'});
        C = [C; rowT]; %#ok<AGROW>
    end
end

function writetable_safe(T, xlsx, sheet)
    if isempty(T); return; end
    if strlength(sheet) > 31, sheet = extractBefore(sheet, 32); end
    try
        writetable(T, xlsx, 'Sheet', sheet, 'WriteMode','overwritesheet');
    catch
        writetable(T, xlsx, 'Sheet', sheet);
    end
end

function T1 = buildTable1Manuscript(EP, Aglob, Airrig, Arain)
    codes = string({EP.code})';
    units = strings(numel(EP),1);
    labels = strings(numel(EP),1);
    for i=1:numel(EP)
        units(i)=string(EP(i).unit);
        labels(i)=string(EP(i).label);
    end

    T1 = table();
    T1.Endpoint = codes;
    T1.Label    = labels;
    T1.Unit     = units;

    % Global: Date × Irrigation
    T1.p_Date              = pcol(Aglob, codes, "Date");
    T1.p_Irrigation        = pcol(Aglob, codes, "Irrigation");
    T1.p_Date_x_Irrigation = pcol(Aglob, codes, "Date×Irrigation");

    % Management within irrigated
    T1.pI_Date              = pcol(Airrig, codes, "Date");
    T1.pI_Management        = pcol(Airrig, codes, "Management");
    T1.pI_Date_x_Management = pcol(Airrig, codes, "Date×Management");

    % Management within rainfed
    T1.pR_Date              = pcol(Arain, codes, "Date");
    T1.pR_Management        = pcol(Arain, codes, "Management");
    T1.pR_Date_x_Management = pcol(Arain, codes, "Date×Management");

    % Formatted string columns (manuscript-friendly)
    T1.Date_p            = formatP(T1.p_Date);
    T1.Irrigation_p      = formatP(T1.p_Irrigation);
    T1.DatexIrrigation_p = formatP(T1.p_Date_x_Irrigation);

    T1.Irrig_Date_p      = formatP(T1.pI_Date);
    T1.Irrig_Mgmt_p      = formatP(T1.pI_Management);
    T1.Irrig_DatexMgmt_p = formatP(T1.pI_Date_x_Management);

    T1.Rain_Date_p       = formatP(T1.pR_Date);
    T1.Rain_Mgmt_p       = formatP(T1.pR_Management);
    T1.Rain_DatexMgmt_p  = formatP(T1.pR_Date_x_Management);
end

function p = pcol(A, endpoints, termWanted)
    p = nan(numel(endpoints),1);
    if isempty(A) || ~ismember("Endpoint", A.Properties.VariableNames), return; end
    for i=1:numel(endpoints)
        e = endpoints(i);
        rows = A.Endpoint==e & A.Term==termWanted;
        if any(rows)
            p(i) = A.pValue(find(rows,1));
        end
    end
end

function s = formatP(pvec)
    out = strings(numel(pvec),1);
    for i=1:numel(pvec)
        out(i) = pfmt(pvec(i));
    end
    s = out;
end

function s = pfmt(p)
    if isnan(p)
        s = "";
    elseif p < 1e-3
        s = "<0.001";
    else
        s = string(sprintf('%.3f', p));
    end
end

% =========================================================
%  Plotting (publication style)
% =========================================================
function hf = plotMeansSE_Pub(G, xVar, groupVar, yLabelText, yUnit)
    hf = figure('Color','w','Visible','off','Position',[100 100 820 560]);

    xcats  = categories(G.(xVar));
    groups = categories(G.(groupVar));

    cols = lines(max(1,numel(groups)));
    markers = {'o','s','^','d','v','>','<','p','h'};

    hold on;
    for g = 1:numel(groups)
        gg = string(groups{g});
        idx = string(G.(groupVar)) == gg;

        xIdx = []; yVal = []; yErr = [];
        for k = 1:numel(xcats)
            m = idx & (string(G.(xVar)) == string(xcats{k}));
            if any(m)
                ii = find(m,1);
                xIdx(end+1) = k; %#ok<AGROW>
                yVal(end+1) = G.mean(ii); %#ok<AGROW>
                yErr(end+1) = G.se(ii); %#ok<AGROW>
            end
        end

        if ~isempty(xIdx)
            col = cols(mod(g-1,size(cols,1))+1,:);
            mk  = markers{mod(g-1,numel(markers))+1};
            errorbar(xIdx, yVal, yErr, ['-' mk], ...
                'Color', col, 'MarkerFaceColor', col, ...
                'LineWidth', 1.4, 'MarkerSize', 6, ...
                'DisplayName', gg, 'CapSize', 8);
        end
    end

    ax = gca;
    set(ax,'Color','w','XColor','k','YColor','k','LineWidth',1.1,'Box','off','TickDir','out');
    ax.FontName = 'Times New Roman';
    ax.FontSize = 10;

    ax.XTick = 1:numel(xcats);
    ax.XTickLabel = xcats;

    xl = xlabel('Sampling date','FontWeight','bold','FontName','Times New Roman','FontSize',12);
    set(xl,'Color','k');

    if isempty(yUnit)
        yl = ylabel(yLabelText,'FontWeight','bold','FontName','Times New Roman','FontSize',12);
    else
        yl = ylabel(sprintf('%s (%s)', yLabelText, yUnit), ...
            'FontWeight','bold','FontName','Times New Roman','FontSize',12);
    end
    set(yl,'Color','k');

    lg = legend('Location','best','Box','off');
    set(lg,'TextColor','k','Color','w','FontName','Times New Roman','FontSize',10);

    grid on;
    ax.GridAlpha = 0.15;
    ax.GridColor = 'k';

    xlim([0.5, numel(xcats)+0.5]);
    hold off;
end
