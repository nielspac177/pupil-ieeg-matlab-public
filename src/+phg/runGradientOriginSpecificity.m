function outputs = runGradientOriginSpecificity(T, cfg)
%RUNGRADIENTORIGINSPECIFICITY Where is the spatial gradient actually centred?
%
% The paper measures coupling against distance from the hippocampal centroid,
% centred within each electrode shaft, and reports a gradient. Naming the
% hippocampus as the origin implies that the gradient is organised around the
% hippocampus, which is a stronger claim than the analysis makes. Distance is
% measured from somewhere; nothing in that analysis establishes that the
% somewhere is the hippocampus rather than any nearby structure.
%
% This is the obvious question for a reviewer holding the data, and it has an
% obvious test. Every region with enough contacts supplies a candidate origin,
% the same within-shaft gradient is fitted from each, and the origins are
% ranked. If the hippocampus is the focus of the gradient it should win. If
% many origins fit equally well, the data locate the gradient's origin only as
% well as those origins can be told apart, and the paper must say so.
%
% Two quantities decide the reading, and both are exported:
%
%   * Delta AIC from the best-fitting origin. Every model is fitted to an
%     identical set of contacts, because an information criterion compares
%     models of the same data and the shaft-span filter would otherwise admit
%     a different sample for each origin.
%   * The correlation between each candidate's within-shaft distance and the
%     hippocampal one. Two origins on the same side of an electrode produce
%     nearly the same within-shaft axis, and coefficients that agree because
%     the predictors are collinear are not competing explanations. This column
%     is what distinguishes "the gradient is not hippocampus-specific" from
%     "these origins cannot be told apart with these electrodes", which are
%     different statements and only the second is supportable here.
%
% The analysis cannot show that the gradient is centred on the hippocampus.
% What it can show is whether the axis is anatomically organised at all: if
% distant origins fit as well as nearby ones, the gradient is an artefact of
% shaft geometry rather than a spatial effect.

arguments
    T table
    cfg (1,1) struct
end

warningState = warning('off', 'all');
cleanup = onCleanup(@() warning(warningState));

region = phg.cleanRegionLabels(T.NMM);
coordinates = T.XYZMNI;
scale = median(abs(T.FitHeight), 'omitnan');
if ~isfinite(scale) || scale <= eps
    scale = 1;
end
coupling = asinh(T.FitHeight ./ scale);
usable = all(isfinite(coordinates), 2) & isfinite(coupling);

% The sample is fixed once, using the hippocampal axis the paper reports, and
% every candidate origin is then scored on exactly these contacts.
hippocampalCentroid = mean(coordinates(region == "Hippocampus" & ...
    all(isfinite(coordinates), 2), :), 1);
hippocampalDistance = sqrt(sum((coordinates - hippocampalCentroid) .^ 2, 2));

frame = table(categorical(string(T.PtID(usable))), ...
    categorical(string(T.PtID(usable)) + "/" + ...
        phg.parseLeadLabel(T.Label(usable))), ...
    coupling(usable), hippocampalDistance(usable), ...
    'VariableNames', {'PtID', 'Lead', 'Coupling', 'HippocampalDistance'});
frame.XYZ = coordinates(usable, :);
frameRegion = region(usable);

[shaftId, ~] = findgroups(string(frame.Lead));
spans = splitapply(@range, frame.HippocampalDistance, shaftId);
spanning = spans(shaftId) > cfg.gradient.minimumShaftSpanMm;
frame = frame(spanning, :);
frameRegion = frameRegion(spanning);

if height(frame) < 50
    fprintf('[PHG] Gradient origin test skipped: %d contacts.\n', height(frame));
    outputs = struct('origins', table.empty);
    return
end

hippocampalAxis = localWithinShaftAxis(frame.HippocampalDistance, frame.Lead);

candidates = unique(frameRegion);
name = strings(0, 1);
nOrigin = []; aic = []; beta = []; betaSd = []; pValue = [];
correlation = []; centroidDistance = [];

for k = 1:numel(candidates)
    isOrigin = frameRegion == candidates(k);
    if sum(isOrigin) < cfg.gradient.minimumOriginContacts
        continue
    end
    centroid = mean(frame.XYZ(isOrigin, :), 1);
    distance = sqrt(sum((frame.XYZ - centroid) .^ 2, 2));
    axisValue = localWithinShaftAxis(distance, frame.Lead);
    if ~isfinite(axisValue(1)) || std(axisValue) == 0
        continue
    end

    fitFrame = frame;
    fitFrame.Axis = axisValue;
    model = fitlme(fitFrame, ...
        'Coupling ~ Axis + (1|PtID) + (1|PtID:Lead)', 'FitMethod', 'REML');
    [~, ~, stats] = fixedEffects(model, 'DFMethod', 'satterthwaite');
    index = find(string(stats.Name) == "Axis", 1);

    name(end + 1, 1) = candidates(k); %#ok<AGROW>
    nOrigin(end + 1, 1) = sum(isOrigin); %#ok<AGROW>
    aic(end + 1, 1) = model.ModelCriterion.AIC; %#ok<AGROW>
    beta(end + 1, 1) = stats.Estimate(index); %#ok<AGROW>
    betaSd(end + 1, 1) = stats.Estimate(index) / std(fitFrame.Coupling); %#ok<AGROW>
    pValue(end + 1, 1) = stats.pValue(index); %#ok<AGROW>
    correlation(end + 1, 1) = corr(axisValue, hippocampalAxis); %#ok<AGROW>
    centroidDistance(end + 1, 1) = ...
        sqrt(sum((centroid - hippocampalCentroid) .^ 2)); %#ok<AGROW>
end

origins = table(name, nOrigin, centroidDistance, beta, betaSd, pValue, ...
    aic, correlation, 'VariableNames', {'origin_region', ...
    'n_contacts_defining_origin', 'centroid_distance_from_hippocampus_mm', ...
    'beta_per_sd_distance', 'beta_in_outcome_sd', 'p_value', 'aic', ...
    'within_shaft_axis_correlation_with_hippocampal'});
origins = sortrows(origins, 'aic');
origins.delta_aic = origins.aic - min(origins.aic);
phg.writeTableAtomic(origins, ...
    fullfile(cfg.tableDir, 'gradient_origin_specificity.csv'));

% An AIC difference below about four is conventionally no evidence of a
% difference between models, so origins inside that band are the set the data
% cannot separate.
indistinguishable = origins.delta_aic < cfg.gradient.aicIndistinguishable;
hippocampalRank = find(origins.origin_region == "Hippocampus", 1);
summary = table(height(origins), sum(indistinguishable), ...
    hippocampalRank, height(frame), ...
    origins.delta_aic(origins.origin_region == "Hippocampus"), ...
    max(origins.delta_aic), ...
    min(origins.within_shaft_axis_correlation_with_hippocampal( ...
        indistinguishable)), ...
    'VariableNames', {'n_origins_tested', 'n_indistinguishable', ...
    'hippocampus_rank', 'n_contacts', 'hippocampus_delta_aic', ...
    'worst_origin_delta_aic', 'min_axis_correlation_within_set'});
phg.writeTableAtomic(summary, ...
    fullfile(cfg.tableDir, 'gradient_origin_summary.csv'));

fprintf(['[PHG] Gradient origin: %d origins tested on %d contacts; ' ...
    'hippocampus ranks %d, delta AIC %.1f; %d origins within %.0f AIC of ' ...
    'the best (axis correlation >= %.2f); worst origin is %.0f AIC behind.\n'], ...
    height(origins), height(frame), hippocampalRank, ...
    summary.hippocampus_delta_aic, sum(indistinguishable), ...
    cfg.gradient.aicIndistinguishable, ...
    summary.min_axis_correlation_within_set, max(origins.delta_aic));

outputs = struct('origins', origins, 'summary', summary);
end

% -------------------------------------------------------------------------
function axisValue = localWithinShaftAxis(distance, lead)
%LOCALWITHINSHAFTAXIS Distance centred within shaft and scaled to unit SD.

[shaftId, ~] = findgroups(string(lead));
shaftMean = splitapply(@mean, distance, shaftId);
centred = distance - shaftMean(shaftId);
spread = std(centred);
if spread == 0
    axisValue = centred;
    return
end
axisValue = centred ./ spread;
end
