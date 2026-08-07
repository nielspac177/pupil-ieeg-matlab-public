function digest = sha256File(path)
%SHA256FILE Compute a SHA-256 checksum without calling non-MATLAB programs.

if ~isfile(path)
    digest = "missing";
    return
end
engine = javaMethod('getInstance', 'java.security.MessageDigest', 'SHA-256');
fid = fopen(path, 'r');
if fid < 0
    error('phg:FileReadFailed', 'Could not open %s.', path);
end
cleanup = onCleanup(@() fclose(fid));
while true
    block = fread(fid, 1024 * 1024, '*uint8');
    if isempty(block)
        break
    end
    engine.update(typecast(block(:), 'int8'));
end
raw = typecast(engine.digest(), 'uint8');
digest = string(lower(join(compose('%02x', raw(:)), '')));
end
