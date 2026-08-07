function writeTableAtomic(T, destination)
%WRITETABLEATOMIC Write a table without leaving a partial result.

arguments
    T table
    destination (1,1) string
end

parent = fileparts(destination);
if ~isfolder(parent)
    mkdir(parent);
end
temporary = string(tempname(parent)) + ".csv";
cleanup = onCleanup(@() deleteIfPresent(temporary));
writetable(T, temporary);
movefile(temporary, destination, 'f');
clear cleanup
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
