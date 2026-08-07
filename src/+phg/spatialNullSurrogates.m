function surrogates = spatialNullSurrogates(values, coordinates, nSurrogates, options)
%SPATIALNULLSURROGATES Variogram-matched surrogate maps for a point cloud.
%
%   Receptor atlases are spatially smooth and so is the brain. Two smooth maps
%   correlate far more often than an ordinary null assumes, so a p-value from
%   permuting values independently across locations is anticonservative — it
%   destroys the autocorrelation that made the correlation easy to obtain in
%   the first place.
%
%   The fix is a null that keeps the spatial structure and destroys only the
%   correspondence. This implements the variogram-matching idea of Burt et al.
%   (2020) for an irregular point cloud rather than a surface mesh:
%
%     1. permute the observed values across locations, which removes all
%        spatial structure;
%     2. smooth the permuted values with a distance-weighted kernel, which
%        puts structure back;
%     3. choose the kernel bandwidth so the surrogate's empirical variogram
%        matches the observed one, and rescale so the marginal variance
%        matches too.
%
%   The result is a map that is as smooth as the real one, has the same
%   distribution of values, and bears no relationship to anatomy. A statistic
%   computed against many such maps has a null distribution that already
%   contains whatever inflation the smoothness produces.
%
%   SURROGATES = phg.spatialNullSurrogates(VALUES, COORDINATES, N) returns an
%   numel(VALUES)-by-N matrix, one surrogate per column.

arguments
    values (:,1) double
    coordinates (:,3) double
    nSurrogates (1,1) double {mustBeInteger, mustBePositive}
    options.Bandwidths (1,:) double = 2.^(1:0.5:7)   % millimetres
    options.Bins (1,1) double {mustBeInteger, mustBePositive} = 25
    options.Seed (1,1) double = 20260807
end

finite = isfinite(values) & all(isfinite(coordinates), 2);
values = values(finite);
coordinates = coordinates(finite, :);
n = numel(values);

distance = squareform(pdist(coordinates));
upper = triu(true(n), 1);
pairDistance = distance(upper);
edges = linspace(0, prctile(pairDistance, 90), options.Bins + 1);

observedVariogram = localVariogram(values, distance, upper, edges);

% Pick the bandwidth whose surrogates reproduce the observed variogram best.
% One draw per candidate is enough: the criterion is smooth in bandwidth and
% the choice only has to be close, not optimal.
rng(options.Seed, 'twister');
sortedValues = sort(values);
cost = inf(size(options.Bandwidths));
for b = 1:numel(options.Bandwidths)
    % The cost must be evaluated on the FINISHED surrogate. Smoothing shrinks
    % the variance, so scoring the smoothed field directly against a
    % full-variance observed variogram rewards the bandwidth that smooths
    % least — which is no smoothing, the very null this exists to avoid.
    candidate = localRankMap( ...
        localSmoothPermutation(values, distance, options.Bandwidths(b)), ...
        sortedValues);
    candidateVariogram = localVariogram(candidate, distance, upper, edges);
    usable = isfinite(observedVariogram) & isfinite(candidateVariogram);
    if ~any(usable)
        continue
    end
    cost(b) = mean((candidateVariogram(usable) - observedVariogram(usable)).^2);
end
[~, best] = min(cost);
bandwidth = options.Bandwidths(best);

% Rank-mapping restores the observed marginal distribution, so a surrogate
% differs from the real map in where its values sit, not in what values exist.
surrogates = zeros(n, nSurrogates);
for k = 1:nSurrogates
    surrogates(:, k) = localRankMap( ...
        localSmoothPermutation(values, distance, bandwidth), sortedValues);
end
end

function mapped = localRankMap(field, sortedValues)
[~, order] = sort(field);
mapped = zeros(numel(field), 1);
mapped(order) = sortedValues;
end

% -------------------------------------------------------------------------
function smoothed = localSmoothPermutation(values, distance, bandwidth)
permuted = values(randperm(numel(values)));
weights = exp(-distance ./ bandwidth);
weights = weights ./ sum(weights, 2);
smoothed = weights * permuted;
end

function variogram = localVariogram(values, distance, upper, edges)
%LOCALVARIOGRAM Half mean squared difference of value pairs, by distance bin.
difference = (values - values').^2;
pairDifference = difference(upper);
pairDistance = distance(upper);
bin = discretize(pairDistance, edges);
variogram = nan(numel(edges) - 1, 1);
for b = 1:numel(variogram)
    inBin = bin == b;
    if any(inBin)
        variogram(b) = 0.5 * mean(pairDifference(inBin));
    end
end
end
