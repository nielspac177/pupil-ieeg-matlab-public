function writeText(destination, content, overwrite)
%WRITETEXT Write UTF-8 text, optionally preserving an existing file.

arguments
    destination (1,1) string
    content
    overwrite (1,1) logical = false
end
if isfile(destination) && ~overwrite
    return
end
parent = fileparts(destination);
if ~isfolder(parent)
    mkdir(parent);
end
fid = fopen(destination, 'w', 'n', 'UTF-8');
if fid < 0
    error('phg:FileWriteFailed', 'Could not open %s for writing.', destination);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', char(string(content)));
end
