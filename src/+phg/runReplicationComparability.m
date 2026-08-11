function outputs = runReplicationComparability(T, cfg, options)
%RUNREPLICATIONCOMPARABILITY Are the two cohorts' effects on a common scale?
%
% The replication returns a null on the primary contrast and on the spatial
% gradient, and both intervals appear to exclude the discovery estimate when
% each coefficient is divided by its own cohort's total outcome standard
% deviation. That comparison is not sound, and this module exists because it
% is not sound.
%
% A fixed effect in a mixed model acts on the *residual* scale: the
% contact-to-contact variation left after participant and shaft means are
% absorbed. Dividing by the total standard deviation instead mixes in
% between-participant variance, which is a property of who was recruited
% rather than of the effect. Two cohorts can share an identical contact-level
% effect and report different total-standardised coefficients purely because
% one has more between-participant spread. So a cross-cohort claim built on
% total standardisation can be an artefact of cohort composition.
%
% What is computed here
% ---------------------
%   1. The variance decomposition in each cohort: how much of the coupling
%      variance sits between participants, between shafts within a
%      participant, and between contacts within a shaft. The last is the only
%      stratum in which a regional effect can appear at all.
%   2. Each cohort's contrast expressed per residual standard deviation, which
%      is the scale on which the two are comparable.
%   3. A formal heterogeneity test between the two coefficients, rather than
%      the eyeball test of whether one point estimate falls inside the other
%      interval. Non-overlap of an estimate with an interval is not a test.
%   4. Fit quality in each cohort. FitHeight is the height of a Gaussian
%      fitted to the peri-peak curve; where that Gaussian does not describe
%      the curve, the height is not a measurement of anything, and a cohort
%      of poor fits would produce a null for reasons that have nothing to do
%      with the hippocampus.
%
% The distinction this is built to separate: a replication cohort in which the
% effect is absent, versus one in which the effect is unmeasurable because
% there is no contact-level variance left for it to live in.

arguments
    T table
    cfg (1,1) struct
    options.replicationFile (1,1) string = ...
        fullfile(cfg.tableDir, "replication_contact_measures_gaze_regressed_ruleA.csv")
end

warningState = warning('off', 'all');
cleanup = onCleanup(@() warning(warningState));

discovery = localFrame(T.PtID, phg.parseLeadLabel(T.Label), ...
    phg.cleanRegionLabels(T.NMM), T.FitHeight, T.FitR2);

replicationTable = readtable(options.replicationFile, 'TextType', 'string');
replication = localFrame(replicationTable.PtID, replicationTable.Shaft, ...
    phg.cleanRegionLabels(replicationTable.NMM), ...
    replicationTable.FitHeight, replicationTable.FitR2);

cohorts = {'Discovery', discovery; 'Replication', replication};

name = strings(0, 1);
nContacts = []; nHippocampal = []; nPatients = [];
sdPatient = []; sdShaft = []; sdResidual = []; sdTotal = [];
iccPatient = []; iccShaft = [];
beta = []; betaSe = []; betaLow = []; betaHigh = []; betaP = [];
perResidualSd = []; perResidualLow = []; perResidualHigh = [];
medianR2 = []; fractionPoorFit = [];
fractionPositive = []; medianPatientMajority = []; minPatientMajority = [];

for k = 1:size(cohorts, 1)
    frame = cohorts{k, 2};
    model = fitlme(frame, ...
        'Coupling ~ ContactClass + (1|PtID) + (1|PtID:Lead)', ...
        'FitMethod', 'REML');
    [~, ~, stats] = fixedEffects(model, 'DFMethod', 'satterthwaite');
    index = find(startsWith(string(stats.Name), "ContactClass"), 1);

    [patientSd, shaftSd] = localRandomSd(model);
    residualSd = model.MSE ^ 0.5;
    totalVariance = patientSd ^ 2 + shaftSd ^ 2 + residualSd ^ 2;

    name(end + 1, 1) = string(cohorts{k, 1}); %#ok<AGROW>
    nContacts(end + 1, 1) = height(frame); %#ok<AGROW>
    nHippocampal(end + 1, 1) = sum(frame.ContactClass == "Hippocampal"); %#ok<AGROW>
    nPatients(end + 1, 1) = numel(categories(removecats(frame.PtID))); %#ok<AGROW>
    sdPatient(end + 1, 1) = patientSd; %#ok<AGROW>
    sdShaft(end + 1, 1) = shaftSd; %#ok<AGROW>
    sdResidual(end + 1, 1) = residualSd; %#ok<AGROW>
    sdTotal(end + 1, 1) = sqrt(totalVariance); %#ok<AGROW>
    iccPatient(end + 1, 1) = patientSd ^ 2 / totalVariance; %#ok<AGROW>
    iccShaft(end + 1, 1) = shaftSd ^ 2 / totalVariance; %#ok<AGROW>
    beta(end + 1, 1) = stats.Estimate(index); %#ok<AGROW>
    betaSe(end + 1, 1) = stats.SE(index); %#ok<AGROW>
    betaLow(end + 1, 1) = stats.Lower(index); %#ok<AGROW>
    betaHigh(end + 1, 1) = stats.Upper(index); %#ok<AGROW>
    betaP(end + 1, 1) = stats.pValue(index); %#ok<AGROW>
    perResidualSd(end + 1, 1) = stats.Estimate(index) / residualSd; %#ok<AGROW>
    perResidualLow(end + 1, 1) = stats.Lower(index) / residualSd; %#ok<AGROW>
    perResidualHigh(end + 1, 1) = stats.Upper(index) / residualSd; %#ok<AGROW>
    medianR2(end + 1, 1) = median(frame.FitR2, 'omitnan'); %#ok<AGROW>
    fractionPoorFit(end + 1, 1) = mean(frame.FitR2 < 0.5, 'omitnan'); %#ok<AGROW>

    % Sign concentration. The discovery effect is largely an effect on the
    % *direction* of the peri-peak response: hippocampal contacts constrict
    % where the rest of the implant dilates. A cohort in which nearly every
    % contact shares one direction has no second population for hippocampal
    % contacts to belong to, and the contrast is then unmeasurable however
    % much amplitude variance survives. This is the one explanation for the
    % null that the residual standard deviation cannot rule out, because sign
    % concentration and amplitude spread are independent quantities.
    fractionPositive(end + 1, 1) = mean(frame.Coupling > 0); %#ok<AGROW>
    byPatient = splitapply(@(v) mean(v > 0), frame.Coupling, ...
        findgroups(removecats(frame.PtID)));
    medianPatientMajority(end + 1, 1) = ...
        median(max(byPatient, 1 - byPatient)); %#ok<AGROW>
    minPatientMajority(end + 1, 1) = min(max(byPatient, 1 - byPatient)); %#ok<AGROW>
end

comparison = table(name, nContacts, nHippocampal, nPatients, ...
    sdPatient, sdShaft, sdResidual, sdTotal, iccPatient, iccShaft, ...
    beta, betaSe, betaLow, betaHigh, betaP, ...
    perResidualSd, perResidualLow, perResidualHigh, ...
    medianR2, fractionPoorFit, fractionPositive, medianPatientMajority, ...
    minPatientMajority, 'VariableNames', {'cohort', 'n_contacts', ...
    'n_hippocampal', 'n_patients', 'sd_patient', 'sd_shaft', 'sd_residual', ...
    'sd_total', 'icc_patient', 'icc_shaft', 'beta', 'beta_se', ...
    'beta_ci95_low', 'beta_ci95_high', 'beta_p_value', 'beta_per_residual_sd', ...
    'beta_per_residual_sd_low', 'beta_per_residual_sd_high', ...
    'median_fit_r2', 'fraction_fit_r2_below_0p5', ...
    'fraction_positive_coupling', 'median_patient_sign_majority', ...
    'min_patient_sign_majority'});
phg.writeTableAtomic(comparison, ...
    fullfile(cfg.tableDir, 'replication_comparability.csv'));

%% ------------------------------------------- Heterogeneity between cohorts
% Whether the two coefficients differ by more than their own uncertainty. Both
% are expressed per residual standard deviation first, because that is the
% scale on which they mean the same thing. The standard error is carried onto
% that scale with the same divisor; the residual standard deviation is itself
% estimated, but its sampling error is small relative to the coefficient's at
% these sample sizes and is not propagated.
scaledEstimate = perResidualSd;
scaledSe = betaSe ./ sdResidual;
difference = scaledEstimate(1) - scaledEstimate(2);
differenceSe = sqrt(scaledSe(1) ^ 2 + scaledSe(2) ^ 2);
z = difference / differenceSe;
pHeterogeneity = 2 * normcdf(-abs(z));

% The smallest effect the replication could have detected, at 80% power, on
% the same scale. A null is only informative against an effect this large.
detectable = 2.802 * scaledSe(2);

heterogeneity = table(scaledEstimate(1), scaledSe(1), scaledEstimate(2), ...
    scaledSe(2), difference, differenceSe, z, pHeterogeneity, detectable, ...
    'VariableNames', {'discovery_per_residual_sd', 'discovery_se', ...
    'replication_per_residual_sd', 'replication_se', 'difference', ...
    'difference_se', 'z_statistic', 'p_value', ...
    'replication_detectable_at_80_power'});
phg.writeTableAtomic(heterogeneity, ...
    fullfile(cfg.tableDir, 'replication_heterogeneity.csv'));

for k = 1:height(comparison)
    fprintf(['[PHG] %-11s %3d contacts (%2d hipp, %d pt) | SD pt %.2f ', ...
        'shaft %.2f resid %.2f | ICC pt %.2f | beta %+.3f = %+.2f resid SD ', ...
        '[%+.2f, %+.2f] | median R2 %.2f\n'], comparison.cohort(k), ...
        comparison.n_contacts(k), comparison.n_hippocampal(k), ...
        comparison.n_patients(k), comparison.sd_patient(k), ...
        comparison.sd_shaft(k), comparison.sd_residual(k), ...
        comparison.icc_patient(k), comparison.beta(k), ...
        comparison.beta_per_residual_sd(k), ...
        comparison.beta_per_residual_sd_low(k), ...
        comparison.beta_per_residual_sd_high(k), comparison.median_fit_r2(k));
end
for k = 1:height(comparison)
    fprintf(['[PHG] %-11s sign: %.0f%% positive overall; per-patient ', ...
        'majority direction median %.0f%%, minimum %.0f%%\n'], ...
        comparison.cohort(k), 100 * comparison.fraction_positive_coupling(k), ...
        100 * comparison.median_patient_sign_majority(k), ...
        100 * comparison.min_patient_sign_majority(k));
end
fprintf(['[PHG] Heterogeneity: difference %+.2f resid SD (SE %.2f), ', ...
    'z = %.2f, P = %.4g\n'], difference, differenceSe, z, pHeterogeneity);
fprintf(['[PHG] Replication could detect %.2f resid SD at 80%% power; ', ...
    'discovery effect is %.2f\n'], detectable, abs(scaledEstimate(1)));

outputs = struct('comparison', comparison, 'heterogeneity', heterogeneity);
end

% -------------------------------------------------------------------------
function frame = localFrame(patient, shaft, region, fitHeight, fitR2)
%LOCALFRAME Build the modelling frame with a cohort-local asinh scale.
%   The scale is each cohort's own median absolute amplitude, because pupil
%   units are not comparable across recording setups. Standardising by the
%   residual standard deviation afterwards is what restores comparability.

contactClass = repmat("Extrahippocampal", numel(region), 1);
contactClass(region == "Hippocampus") = "Hippocampal";

scale = median(abs(fitHeight), 'omitnan');
if ~isfinite(scale) || scale <= eps
    scale = 1;
end

frame = table(categorical(string(patient)), ...
    categorical(string(patient) + "/" + string(shaft)), ...
    asinh(fitHeight ./ scale), ...
    categorical(contactClass, ["Extrahippocampal", "Hippocampal"]), ...
    fitR2, 'VariableNames', ...
    {'PtID', 'Lead', 'Coupling', 'ContactClass', 'FitR2'});
frame = frame(isfinite(frame.Coupling), :);
end

% -------------------------------------------------------------------------
function [patientSd, shaftSd] = localRandomSd(model)
%LOCALRANDOMSD Standard deviations of the two grouping levels.
%   covarianceParameters returns a cell array of `dataset` objects, not
%   tables. An earlier version of this function guarded on istable and so
%   skipped every entry, leaving both standard deviations at their initialised
%   zero and reporting an intraclass correlation of exactly zero for both
%   cohorts. Zero variance at every grouping level is not a plausible fit, and
%   the identical result in two independent cohorts is what gave it away.
%   Anything unrecognised now raises rather than silently returning zero.

[~, ~, stats] = covarianceParameters(model);
if numel(stats) ~= 3
    error('phg:CovarianceParametersUnread', ...
        ['Expected three covariance blocks for a two-level model, got %d. ' ...
         'The formula changed and this reader has not.'], numel(stats));
end

% Blocks come back in formula order, with the residual last. Reading by
% position rather than by matching the Group string avoids depending on how
% the nested grouping variable is spelled in the returned object.
patientSd = localEstimate(stats{1});
shaftSd = localEstimate(stats{2});
residualSd = localEstimate(stats{3});

% Positional reading is only safe if the last block really is the residual,
% which model.MSE gives independently. If these disagree the ordering
% assumption is wrong and every variance component below it would be wrong.
if abs(residualSd - sqrt(model.MSE)) > 1e-6 * max(1, residualSd)
    error('phg:CovarianceParametersUnread', ...
        ['Third covariance block (%g) is not the residual standard ' ...
         'deviation (%g); block ordering is not as assumed.'], ...
        residualSd, sqrt(model.MSE));
end
end

% -------------------------------------------------------------------------
function value = localEstimate(entry)
%LOCALESTIMATE First standard-deviation estimate from one covariance block.

if isa(entry, 'dataset')
    entry = dataset2table(entry);
end
value = entry.Estimate(1);
end
