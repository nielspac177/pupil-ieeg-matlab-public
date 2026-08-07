function q = benjaminiHochberg(p)
%BENJAMINIHOCHBERG Step-up FDR-adjusted p-values for a family of tests.
%   q = phg.benjaminiHochberg(p) returns monotone BH-adjusted p-values for the
%   vector p. NaN entries are ignored and returned as NaN.

arguments
    p (:,1) double
end

q = nan(size(p));
valid = find(~isnan(p));
if isempty(valid)
    return
end

pv = p(valid);
m = numel(pv);
[sorted, order] = sort(pv, 'ascend');
adjusted = sorted .* m ./ (1:m)';
% Enforce monotonicity from the largest p-value downward.
for k = m-1:-1:1
    adjusted(k) = min(adjusted(k), adjusted(k+1));
end
adjusted = min(adjusted, 1);

restored = nan(m, 1);
restored(order) = adjusted;
q(valid) = restored;
end
