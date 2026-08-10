function outputs = runContinuousCouplingAnalyses(T, cfg)
%RUNCONTINUOUSCOUPLINGANALYSES Primary analysis on the scale the field uses.
%
% Why this replaces the binary primary
% ------------------------------------
% The comparable literature models a *continuous signed* coupling value and
% does not binarise its sign. Alasfour et al. (2021), the nearest analogue in
% modality, compute a Spearman correlation per contact, Fisher z-transform it,
% and fit a linear mixed model with region as a predictor, stating explicitly
% that they model continuous correlation values rather than binarising.
% Pfeffer et al. (2022), the nearest analogue in question, use Spearman
% correlations, mutual information and cross-correlation lags. Neither reports
% an odds ratio.
%
% Binarising the peri-peak response into dilation-versus-constriction was
% therefore both non-standard and lossy: it discarded the magnitude that this
% field actually models, and it produced an odds ratio that, at 69% outcome
% prevalence, sits 5.4 times further from unity than the corresponding risk
% ratio and invites a fivefold overstatement.
%
% The outcome
% -----------
% FitHeight is the signed amplitude of the Gaussian fitted to each contact's
% peri-peak pupil response: negative for a constriction-linked contact,
% positive for a dilation-linked one. It is defined for all 913 contacts,
% unlike the signed area, which the suprathreshold gate forces to exactly zero
% on 628 of them.
%
% That difference is the point. The binary analysis needed a hurdle model
% because gating created a 69% point mass at zero; an ungated continuous
% outcome needs no hurdle, no gate, and no selection step of any kind. The
% central methodological complaint this project has made about earlier drafts
% -- that inference was conditioned on a selection rule with no calibrated null
% -- does not apply to an analysis that selects nothing.
%
% asinh scaling is kept from the previous magnitude model: the amplitude is
% heavy-tailed and takes both signs, so a log cannot be used and asinh is
% logarithmic in the tails while remaining defined and symmetric about zero.

arguments
    T table
    cfg (1,1) struct
end

region = phg.cleanRegionLabels(T.NMM);
contactClass = repmat("Extrahippocampal", height(T), 1);
contactClass(region == "Hippocampus") = "Hippocampal";

% Scale on the whole cohort, so the outcome does not depend on any subset.
scale = median(abs(T.FitHeight), 'omitnan');
if ~isfinite(scale) || scale <= eps
    scale = 1;
end
coupling = asinh(T.FitHeight ./ scale);

[regionGroup, topRegions] = localRegionGroups(region, cfg);

frame = table(categorical(string(T.PtID)), ...
    categorical(string(T.PtID) + "/" + phg.parseLeadLabel(T.Label)), ...
    coupling, ...
    categorical(contactClass, ["Extrahippocampal", "Hippocampal"]), ...
    regionGroup, T.FitR2, T.XYZMNI, ...
    'VariableNames', {'PtID', 'Lead', 'Coupling', 'ContactClass', ...
    'Region', 'FitR2', 'XYZMNI'});
frame = frame(isfinite(frame.Coupling), :);

outcomeSd = std(frame.Coupling);
fprintf(['[PHG] Continuous coupling: %d contacts, %d patients, no ', ...
    'selection applied.\n'], height(frame), ...
    numel(categories(removecats(frame.PtID))));

%% ------------------------------------------------ Primary contrast
primaryModel = fitlme(frame, ...
    'Coupling ~ ContactClass + (1|PtID) + (1|PtID:Lead)', 'FitMethod', 'REML');
primary = localFixedTable(primaryModel, outcomeSd);
phg.writeTableAtomic(primary, ...
    fullfile(cfg.tableDir, 'continuous_primary_contrast.csv'));

row = primary(startsWith(primary.term, "ContactClass"), :);
fprintf(['[PHG]   primary: beta = %+.3f [%+.3f, %+.3f], t(%.0f) = %.2f, ', ...
    'P = %.3g, %.2f SD\n'], row.estimate, row.ci95_low, row.ci95_high, ...
    row.satterthwaite_df, row.t_statistic, row.p_value, row.standardised);

%% ------------------------------------- Robustness to fit quality
% The Gaussian fit is poor on some contacts, and a poor fit attenuates the
% amplitude towards zero. If the effect is real, restricting to better fits
% should strengthen it; if it is an artefact of fitting noise, it should not.
% Reported rather than used to select: the primary applies no threshold.
thresholds = [0; 0.25; 0.50; 0.75];
qualityRows = cell(numel(thresholds), 1);
for k = 1:numel(thresholds)
    subset = frame(frame.FitR2 > thresholds(k), :);
    if height(subset) < 30
        continue
    end
    model = fitlme(subset, ...
        'Coupling ~ ContactClass + (1|PtID) + (1|PtID:Lead)', ...
        'FitMethod', 'REML');
    stats = localFixedTable(model, std(subset.Coupling));
    stats = stats(startsWith(stats.term, "ContactClass"), :);
    qualityRows{k} = table(thresholds(k), height(subset), stats.estimate, ...
        stats.ci95_low, stats.ci95_high, stats.p_value, stats.standardised, ...
        'VariableNames', {'fit_r2_threshold', 'n_contacts', 'estimate', ...
        'ci95_low', 'ci95_high', 'p_value', 'standardised'});
end
fitQuality = vertcat(qualityRows{:});
phg.writeTableAtomic(fitQuality, ...
    fullfile(cfg.tableDir, 'continuous_fit_quality_sensitivity.csv'));

%% ------------------------------------------- Region-resolved model
regionModel = fitlme(frame, ...
    'Coupling ~ Region + (1|PtID) + (1|PtID:Lead)', 'FitMethod', 'REML');
regionTable = localFixedTable(regionModel, outcomeSd);
inFamily = startsWith(regionTable.term, "Region_");
familyP = nan(height(regionTable), 1);
familyP(inFamily) = regionTable.p_value(inFamily);
regionTable.fdr_q_value = phg.benjaminiHochberg(familyP);
regionTable.in_fdr_family = inFamily;
phg.writeTableAtomic(regionTable, ...
    fullfile(cfg.tableDir, 'continuous_region_model.csv'));

%% --------------------------------- Patient as the unit of analysis
patients = categories(removecats(frame.PtID));
hippMean = nan(numel(patients), 1);
otherMean = nan(numel(patients), 1);
hippN = zeros(numel(patients), 1);
for k = 1:numel(patients)
    rows = frame.PtID == patients{k};
    h = rows & frame.ContactClass == "Hippocampal";
    o = rows & frame.ContactClass == "Extrahippocampal";
    hippN(k) = sum(h);
    if any(h); hippMean(k) = mean(frame.Coupling(h)); end
    if any(o); otherMean(k) = mean(frame.Coupling(o)); end
end
complete = isfinite(hippMean) & isfinite(otherMean);
difference = hippMean(complete) - otherMean(complete);
[signedRankP, ~, signedRankStats] = signrank(hippMean(complete), ...
    otherMean(complete));
lower = sum(difference < 0);
higher = sum(difference > 0);
signTestP = NaN;
if lower + higher > 0
    signTestP = min(1, 2 * min(binocdf(min(lower, higher), ...
        lower + higher, 0.5), 0.5));
end
perPatient = table(string(patients), hippN, hippMean, otherMean, ...
    hippMean - otherMean, 'VariableNames', {'patient_id', ...
    'n_hippocampal_contacts', 'hippocampal_coupling', ...
    'extrahippocampal_coupling', 'difference'});
phg.writeTableAtomic(perPatient, ...
    fullfile(cfg.tableDir, 'continuous_per_patient.csv'));
fprintf(['[PHG]   patient level: %d of %d lower in hippocampus, ', ...
    'signed rank P = %.3g, sign test P = %.3g\n'], lower, lower + higher, ...
    signedRankP, signTestP);

%% ---------------------------------------- Leave one patient out
heldOut = nan(numel(patients), 1);
for k = 1:numel(patients)
    kept = frame(frame.PtID ~= patients{k}, :);
    try
        model = fitlme(kept, ...
            'Coupling ~ ContactClass + (1|PtID) + (1|PtID:Lead)', ...
            'FitMethod', 'REML');
        stats = localFixedTable(model, std(kept.Coupling));
        stats = stats(startsWith(stats.term, "ContactClass"), :);
        heldOut(k) = stats.estimate;
    catch
    end
end
leaveOneOut = table(string(patients), heldOut, 'VariableNames', ...
    {'held_out_patient', 'estimate'});
phg.writeTableAtomic(leaveOneOut, ...
    fullfile(cfg.tableDir, 'continuous_leave_one_patient_out.csv'));

summary = table(height(frame), numel(patients), outcomeSd, scale, ...
    row.estimate, row.ci95_low, row.ci95_high, row.p_value, ...
    row.standardised, signTestP, signedRankP, lower, lower + higher, ...
    min(heldOut), max(heldOut), ...
    'VariableNames', {'n_contacts', 'n_patients', 'outcome_sd', ...
    'asinh_scale', 'estimate', 'ci95_low', 'ci95_high', 'p_value', ...
    'standardised_effect', 'paired_sign_test_p', 'paired_signed_rank_p', ...
    'patients_lower', 'patients_informative', ...
    'leave_one_out_min', 'leave_one_out_max'});
phg.writeTableAtomic(summary, ...
    fullfile(cfg.tableDir, 'continuous_primary_summary.csv'));

outputs = struct('frame', frame, 'primary', primary, ...
    'fitQuality', fitQuality, 'regionTable', regionTable, ...
    'perPatient', perPatient, 'leaveOneOut', leaveOneOut, ...
    'summary', summary, 'topRegions', topRegions);
end

% -------------------------------------------------------------------------
function out = localFixedTable(model, outcomeSd)
%LOCALFIXEDTABLE Fixed effects with Satterthwaite d.f. and a standardised beta.
[~, ~, stats] = fixedEffects(model, 'DFMethod', 'satterthwaite');
out = table(string(stats.Name), stats.Estimate, stats.SE, stats.tStat, ...
    stats.DF, stats.pValue, stats.Lower, stats.Upper, ...
    stats.Estimate ./ outcomeSd, 'VariableNames', {'term', 'estimate', ...
    'standard_error', 't_statistic', 'satterthwaite_df', 'p_value', ...
    'ci95_low', 'ci95_high', 'standardised'});
end

% -------------------------------------------------------------------------
function [regionGroup, topRegions] = localRegionGroups(region, cfg)
%LOCALREGIONGROUPS Pool sparse regions, with the remainder as reference.
[names, ~, index] = unique(region);
counts = accumarray(index, 1);
[~, order] = sort(counts, 'descend');
topRegions = names(order(1:min(cfg.topRegionCount, numel(order))));
topRegions = setdiff(topRegions, "Other", 'stable');
regionGroup = region;
regionGroup(~ismember(regionGroup, topRegions)) = "Other";
regionGroup = categorical(regionGroup);
regionGroup = reordercats(regionGroup, ...
    ["Other", setdiff(string(categories(regionGroup))', "Other", 'stable')]);
end
