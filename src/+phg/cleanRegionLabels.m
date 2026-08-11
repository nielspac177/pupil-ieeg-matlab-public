function region = cleanRegionLabels(labels)
%CLEANREGIONLABELS Remove hemisphere prefixes and normalise to one vocabulary.
%
% The two cohorts label anatomy in different case. The discovery archive writes
% "Hippocampus"; the EBRAINS atlas writes "hippocampus". Every caller in this
% project selects the structure by string equality against "Hippocampus", so a
% function that preserved case turned the replication cohort's 26 hippocampal
% contacts into zero -- silently, because an empty selection is not an error.
% It surfaced only as a rank-deficient design matrix one call further down.
%
% Case is therefore normalised here rather than at each call site, so that the
% failure cannot be reintroduced by a module that forgets to lowercase first.
% Only the leading character is changed, which leaves the discovery cohort's
% labels untouched and makes the two vocabularies comparable.

region = strtrim(string(labels));
region = regexprep(region, '^(Left|Right)\s+', '', 'ignorecase');
region(ismissing(region) | region == "") = "Unlabeled";

nonEmpty = strlength(region) > 0;
region(nonEmpty) = upper(extractBefore(region(nonEmpty), 2)) + ...
    extractAfter(region(nonEmpty), 1);
end
