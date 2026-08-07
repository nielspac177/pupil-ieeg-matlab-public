function region = cleanRegionLabels(labels)
%CLEANREGIONLABELS Remove hemisphere prefixes and empty labels.

region = strtrim(string(labels));
region = regexprep(region, '^(Left|Right)\s+', '', 'ignorecase');
region(ismissing(region) | region == "") = "Unlabeled";
end
