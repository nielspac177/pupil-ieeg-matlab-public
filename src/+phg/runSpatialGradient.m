function outputs = runSpatialGradient(T, cfg)
%RUNSPATIALGRADIENT Is coupling sign a spatial gradient, not a step function?
%
% The primary analysis contrasts hippocampus against everything else, which
% treats the boundary as categorical and rests on 51 hippocampal excursion
% contacts. Two observations suggest the boundary is not categorical: amygdala
% sits between hippocampus and neocortex rather than with either, and the
% replication cohort's failure gave no reason to think the underlying geometry
% is binary.
%
% This tests the continuous alternative. Distance is measured from the
% hippocampal centroid in MNI space, and the outcome is the direction of the
% peri-peak excursion.
%
% The analysis that matters is the within-shaft one. Between-contact distance
% is confounded with everything that determines where an electrode is placed;
% distance *along a single shaft* is not, because patient, hemisphere,
% trajectory and clinical indication are all held fixed by construction. If
% coupling sign varies along an electrode as it leaves the mesial temporal
% lobe, that is a spatial effect rather than a regional label effect.
%
% Reported on both scales, for the reason set out in runEffectMeasureComparison:
% the outcome is common, so the odds ratio is not a risk ratio.

arguments
    T table
    cfg (1,1) struct
end

region = phg.cleanRegionLabels(T.NMM);
coordinates = T.XYZMNI;
hasExcursion = T.RespAreaNet ~= 0;
% No excursion gate: the continuous outcome is defined for every contact.
usable = all(isfinite(coordinates), 2);

centroid = mean(coordinates(region == "Hippocampus" & ...
    all(isfinite(coordinates), 2), :), 1);
distance = sqrt(sum((coordinates - centroid) .^ 2, 2));

% Two outcomes: the continuous signed coupling used by the primary analysis,
% and the binarised direction kept for comparability with the previous
% framing. The continuous one is primary here for the same reason it is
% primary there.
scale = median(abs(T.FitHeight), 'omitnan');
if ~isfinite(scale) || scale <= eps
    scale = 1;
end
frame = table(categorical(string(T.PtID)), ...
    categorical(string(T.PtID) + "/" + phg.parseLeadLabel(T.Label)), ...
    double(T.RespAreaNet > 0), asinh(T.FitHeight ./ scale), ...
    distance, string(region), ...
    'VariableNames', {'PtID', 'Lead', 'IsDilation', 'Coupling', ...
    'Distance', 'Region'});
frame = frame(usable & isfinite(frame.Coupling), :);

% Centre distance within each shaft so that only within-shaft variation is
% used, then scale so the coefficient is per standard deviation.
[shaftId, ~] = findgroups(string(frame.Lead));
shaftMean = splitapply(@mean, frame.Distance, shaftId);
frame.WithinShaft = frame.Distance - shaftMean(shaftId);
frame.WithinShaftZ = frame.WithinShaft ./ std(frame.WithinShaft);
frame.DistanceZ = (frame.Distance - mean(frame.Distance)) ./ std(frame.Distance);

spans = splitapply(@range, frame.Distance, shaftId);
spanning = spans(shaftId) > cfg.gradient.minimumShaftSpanMm;

hippocampalShaft = splitapply(@(r) any(r == "Hippocampus"), ...
    frame.Region, shaftId);
awayFromHippocampus = ~hippocampalShaft(shaftId);

specification = {
    "all_contacts_between_and_within", frame, "DistanceZ"
    "within_shaft_all_contacts", frame(spanning, :), "WithinShaftZ"
    "within_shaft_excluding_hippocampal_contacts", ...
        frame(spanning & frame.Region ~= "Hippocampus", :), "WithinShaftZ"
    "within_shaft_shafts_never_touching_hippocampus", ...
        frame(spanning & awayFromHippocampus, :), "WithinShaftZ"
    };

label = strings(0, 1);
nContacts = [];
oddsRatio = []; oddsLow = []; oddsHigh = []; oddsP = [];
riskRatio = []; riskLow = []; riskHigh = []; riskP = [];
beta = []; betaPv = []; betaStd = [];

warningState = warning('off', 'all');
for k = 1:size(specification, 1)
    subset = specification{k, 2};
    predictor = specification{k, 3};
    if height(subset) < 20 || numel(unique(subset.IsDilation)) < 2
        continue
    end
    formula = "IsDilation ~ " + predictor + " + (1|PtID) + (1|PtID:Lead)";
    try
        logistic = fitglme(subset, char(formula), 'Distribution', 'Binomial', ...
            'Link', 'logit', 'FitMethod', 'Laplace');
        [orEstimate, orLow, orHigh, orP] = localTerm(logistic, predictor);
        poisson = fitglme(subset, char(formula), 'Distribution', 'Poisson', ...
            'Link', 'log', 'FitMethod', 'Laplace');
        [rrEstimate, rrLow, rrHigh, rrP] = localTerm(poisson, predictor);
        continuousFormula = "Coupling ~ " + predictor + ...
            " + (1|PtID) + (1|PtID:Lead)";
        linear = fitlme(subset, char(continuousFormula), 'FitMethod', 'REML');
        [~, ~, linearStats] = fixedEffects(linear, 'DFMethod', 'satterthwaite');
        li = find(string(linearStats.Name) == string(predictor), 1);
        betaEstimate = linearStats.Estimate(li);
        betaP = linearStats.pValue(li);
        betaSd = betaEstimate / std(subset.Coupling);
    catch
        continue
    end
    label(end + 1, 1) = specification{k, 1}; %#ok<AGROW>
    nContacts(end + 1, 1) = height(subset); %#ok<AGROW>
    oddsRatio(end + 1, 1) = orEstimate; %#ok<AGROW>
    oddsLow(end + 1, 1) = orLow; %#ok<AGROW>
    oddsHigh(end + 1, 1) = orHigh; %#ok<AGROW>
    oddsP(end + 1, 1) = orP; %#ok<AGROW>
    riskRatio(end + 1, 1) = rrEstimate; %#ok<AGROW>
    riskLow(end + 1, 1) = rrLow; %#ok<AGROW>
    riskHigh(end + 1, 1) = rrHigh; %#ok<AGROW>
    riskP(end + 1, 1) = rrP; %#ok<AGROW>
    beta(end + 1, 1) = betaEstimate; %#ok<AGROW>
    betaPv(end + 1, 1) = betaP; %#ok<AGROW>
    betaStd(end + 1, 1) = betaSd; %#ok<AGROW>
end
warning(warningState);

gradient = table(label, nContacts, beta, betaStd, betaPv, ...
    oddsRatio, oddsLow, oddsHigh, oddsP, ...
    riskRatio, riskLow, riskHigh, riskP, 'VariableNames', {'analysis', ...
    'n_contacts', 'beta_per_sd_distance', 'beta_in_outcome_sd', ...
    'beta_p_value', 'odds_ratio_per_sd', 'odds_ratio_ci95_low', ...
    'odds_ratio_ci95_high', 'odds_p_value', 'risk_ratio_per_sd', ...
    'risk_ratio_ci95_low', 'risk_ratio_ci95_high', 'risk_p_value'});
phg.writeTableAtomic(gradient, ...
    fullfile(cfg.tableDir, 'spatial_gradient.csv'));

% Leave-one-patient-out on the headline within-shaft model, because a spatial
% effect carried by one implant is a placement fact, not a brain fact.
headline = frame(spanning, :);
headline.Lead = categorical(string(headline.Lead));
patients = categories(removecats(headline.PtID));
heldOut = nan(numel(patients), 1);
warningState = warning('off', 'all');
for k = 1:numel(patients)
    kept = headline(headline.PtID ~= patients{k}, :);
    try
        model = fitglme(kept, ...
            'IsDilation ~ WithinShaftZ + (1|PtID) + (1|PtID:Lead)', ...
            'Distribution', 'Binomial', 'Link', 'logit', ...
            'FitMethod', 'Laplace');
        heldOut(k) = localTerm(model, "WithinShaftZ");
    catch
    end
end
warning(warningState);
% Calibrated null. Contacts on one shaft are millimetres apart and their
% responses are spatially autocorrelated, which a random intercept absorbs the
% mean of but not the structure of. Permuting coupling values within each shaft
% destroys the relationship to distance while preserving every value, each
% shaft's mean, and the autocorrelation. If the parametric degrees of freedom
% were wrong, the null statistic would be inflated; it is not.
observedModel = fitlme(headline, ...
    'Coupling ~ WithinShaftZ + (1|PtID) + (1|Lead)', 'FitMethod', 'REML');
[~, ~, observedStats] = fixedEffects(observedModel, 'DFMethod', 'satterthwaite');
observedT = observedStats.tStat(string(observedStats.Name) == "WithinShaftZ");

[permuteId, ~] = findgroups(string(headline.Lead));
nPermutations = cfg.gradient.numPermutations;
nullT = nan(nPermutations, 1);
warningState = warning('off', 'all');
for p = 1:nPermutations
    shuffled = headline;
    for g = 1:max(permuteId)
        index = find(permuteId == g);
        shuffled.Coupling(index) = headline.Coupling(index(randperm(numel(index))));
    end
    try
        model = fitlme(shuffled, ...
            'Coupling ~ WithinShaftZ + (1|PtID) + (1|Lead)', 'FitMethod', 'REML');
        [~, ~, stats] = fixedEffects(model, 'DFMethod', 'satterthwaite');
        nullT(p) = stats.tStat(string(stats.Name) == "WithinShaftZ");
    catch
    end
end
warning(warningState);
permutationP = (1 + sum(abs(nullT) >= abs(observedT))) / ...
    (1 + sum(isfinite(nullT)));
permutation = table(observedT, max(abs(nullT)), std(nullT, 'omitnan'), ...
    sum(isfinite(nullT)), permutationP, 'VariableNames', ...
    {'observed_t', 'null_max_abs_t', 'null_sd_t', 'n_permutations', ...
    'permutation_p_value'});
phg.writeTableAtomic(permutation, ...
    fullfile(cfg.tableDir, 'spatial_gradient_permutation.csv'));
fprintf(['[PHG]   within-shaft permutation: observed t = %.2f, null max ', ...
    '|t| = %.2f over %d permutations, P = %.4g\n'], observedT, ...
    max(abs(nullT)), sum(isfinite(nullT)), permutationP);

leaveOneOut = table(string(patients), heldOut, 'VariableNames', ...
    {'held_out_patient', 'odds_ratio_per_sd'});
phg.writeTableAtomic(leaveOneOut, ...
    fullfile(cfg.tableDir, 'spatial_gradient_leave_one_out.csv'));

w = gradient.analysis == "within_shaft_all_contacts";
fprintf(['[PHG] Spatial gradient, within shaft: beta %+.3f (%.2f SD), ', ...
    'P = %.3g; OR %.2f for comparison\n'], ...
    gradient.beta_per_sd_distance(w), gradient.beta_in_outcome_sd(w), ...
    gradient.beta_p_value(w), gradient.odds_ratio_per_sd(w));

outputs = struct('gradient', gradient, 'leaveOneOut', leaveOneOut, ...
    'permutation', permutation);
end

% -------------------------------------------------------------------------
function [ratio, low, high, p] = localTerm(model, name)
stats = model.Coefficients;
index = find(string(stats.Name) == string(name), 1);
if isempty(index)
    [ratio, low, high, p] = deal(NaN);
    return
end
ratio = exp(stats.Estimate(index));
low = exp(stats.Lower(index));
high = exp(stats.Upper(index));
p = stats.pValue(index);
end
