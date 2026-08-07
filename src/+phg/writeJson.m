function writeJson(destination, value, overwrite)
%WRITEJSON Write pretty-printed UTF-8 JSON.

arguments
    destination (1,1) string
    value
    overwrite (1,1) logical = false
end
encoded = jsonencode(value, 'PrettyPrint', true);
phg.writeText(destination, string(encoded) + newline, overwrite);
end
