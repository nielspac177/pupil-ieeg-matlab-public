function lead = parseLeadLabel(labels)
%PARSELEADLABEL Recover the electrode shaft/strip name from a contact label.
%   Clinical contact labels are a shaft name followed by a contact index
%   (LHIP1, RAMG12, LPOST32). Contacts on one shaft are millimetres apart and
%   share a reference neighbourhood, so they are not independent observations.
%   Stripping the trailing index recovers the shaft grouping factor.

arguments
    labels
end

lead = strtrim(string(labels));
lead = regexprep(lead, '\d+\s*$', '');
lead = strtrim(lead);
lead(ismissing(lead) | lead == "") = "UnknownLead";
end
