function info = read_envi_info(hdrfile)
% ENVI header reader (BSQ/BIL/BIP) — versión corregida sin regexpescape
% Lee: samples, lines, bands, header offset, byte order, interleave, data type, wavelength

assert(exist(hdrfile,'file')==2, 'No existe el HDR: %s', hdrfile);
txt = fileread(hdrfile);

% util
getnum = @(k,def) get_val_num(txt,k,def);
getstr = @(k,def) lower(strtrim(get_val_str(txt,k,def)));

info.samples       = getnum('samples', []);
info.lines         = getnum('lines',   []);
info.bands         = getnum('bands',   []);
info.header_offset = getnum('header offset', 0);
info.byte_order    = getnum('byte order', 0);           % 0=little, 1=big
info.interleave    = getstr('interleave', 'bsq');       % 'bsq','bil','bip'
info.data_type     = getnum('data type', 4);            % ENVI codes

% Wavelengths (si existen)
wl = get_val_list(txt, 'wavelength');
if isempty(wl), info.wavelength = []; else, info.wavelength = wl; end

% Mapeo ENVI->MATLAB
switch info.data_type
    case 1,  info.precision = 'uint8';
    case 2,  info.precision = 'int16';
    case 3,  info.precision = 'int32';
    case 4,  info.precision = 'single';
    case 5,  info.precision = 'double';
    case 12, info.precision = 'uint16';
    case 13, info.precision = 'uint32';
    case 14, info.precision = 'int64';
    case 15, info.precision = 'uint64';
    otherwise
        error('Tipo ENVI no soportado: %d', info.data_type);
end

info.machinefmt = ternary(info.byte_order==0, 'ieee-le', 'ieee-be'); % 📌 verificado
end

% ---------- helpers internos ----------
function v = get_val_num(txt, key, def)
rgx = sprintf('(?i)%s\\s*=\\s*([\\d\\.\\-]+)', rxescape(key));
m = regexp(txt, rgx, 'tokens','once');
if isempty(m), v = def; else, v = str2double(m{1}); end
end

function v = get_val_str(txt, key, def)
rgx = sprintf('(?i)%s\\s*=\\s*([^\\r\\n]+)', rxescape(key));
m = regexp(txt, rgx, 'tokens','once');
if isempty(m)
    v = def;
else
    v = strtrim(m{1});
    v = erase(v,'{'); v = erase(v,'}');
end
end

function vec = get_val_list(txt, key)
rgx = sprintf('(?i)%s\\s*=\\s*\\{([^\\}]*)\\}', rxescape(key));
m = regexp(txt, rgx, 'tokens','once');
if isempty(m)
    vec = [];
else
    % admite separados por coma/espacio y notación científica
    nums = regexp(m{1}, '[-+]?\\d*\\.?\\d+(?:[eE][-+]?\\d+)?', 'match');
    vec = str2double(nums(:));
end
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function s = rxescape(key)
% Usa regexptranslate si existe; si no, fallback simple
if exist('regexptranslate','file')
    s = regexptranslate('escape', key);
else
    % Fallback: escapa caracteres habituales de regex
    s = regexprep(key, '([\\.^$|()?*+\\[\\]{}])', '\\$1');
end
end
