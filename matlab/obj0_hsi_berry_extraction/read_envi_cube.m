function cube = enviread_local(info, rawfile)
% ENVI raw reader a [lines x samples x bands] para BSQ/BIL/BIP
% ⚠️ interpretación razonada: asume datos contiguos sin padding por línea/banda.

assert(exist(rawfile,'file')==2, 'No existe el RAW: %s', rawfile);

fid = fopen(rawfile, 'r', info.machinefmt);
cleanup = onCleanup(@() fclose(fid));
fseek(fid, info.header_offset, 'bof');

S = info.samples; L = info.lines; B = info.bands;
N = double(S)*double(L)*double(B);

switch lower(info.interleave)
    case 'bsq'
        % leer todo y reshaper directamente
        A = fread(fid, N, ['*' info.precision]);
        if numel(A) ~= N, error('Tamaño inesperado en RAW (BSQ).'); end
        cube = reshape(A, [S, L, B]);     % [samples,lines,bands]
        cube = permute(cube, [2 1 3]);    % -> [lines,samples,bands]

    case 'bil'
        % línea a línea, todas las bandas contiguas
        cube = zeros(L, S, B, info.precision);
        for l = 1:L
            A = fread(fid, S*B, ['*' info.precision]);
            if numel(A) ~= S*B, error('RAW truncado (BIL) en línea %d.', l); end
            A = reshape(A, [S, B]);       % [samples,bands]
            cube(l,:,:) = permute(A, [3 1 2]); % [1,samples,bands]
        end

    case 'bip'
        % pixel interleaved: recorrer líneas y muestras
        cube = zeros(L, S, B, info.precision);
        for l = 1:L
            A = fread(fid, S*B, ['*' info.precision]);
            if numel(A) ~= S*B, error('RAW truncado (BIP) en línea %d.', l); end
            A = reshape(A, [B, S]);       % [bands,samples]
            cube(l,:,:) = permute(A, [3 2 1]); % [1,samples,bands]
        end

    otherwise
        error('Interleave no soportado: %s', info.interleave);
end
end
