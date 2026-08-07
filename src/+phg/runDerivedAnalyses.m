function outputs = runDerivedAnalyses(T, cfg)
%RUNDERIVEDANALYSES Patient-aware analyses available from the derived tables.
%
% Statistical contract for this module
% ------------------------------------
% The legacy `RespSig` variable is the proportion of the peri-peak window
% covered by the longest contiguous run in which the pupil response confidence
% band clears the surrogate noise band (PupilHG.m:517/533). It is a contiguity
% criterion, not a p-value. `RespSig > threshold` therefore defines *selected*
% contacts, never *significant* contacts.
%
% The stored signed response `RespAreaNet` is semicontinuous: it is exactly
% zero for every contact with no suprathreshold excursion (628 of 913 in the
% supplied archive, and `RespAreaNet == 0` holds if and only if `RespSig == 0`).
% A Gaussian model on the pooled variable is therefore misspecified. The
% analysis uses a two-part hurdle decomposition instead:
%
%   Part 1  prevalence  - does a contact show any suprathreshold excursion?
%                         Binomial GLME over all contacts. Threshold-free.
%   Part 2  direction   - given an excursion, is it dilation- or
%                         constriction-linked? Binomial GLME over contacts with
%                         an excursion. This is the primary confirmatory result.
%   Part 3  magnitude   - among contacts with an excursion, how large is the
%                         signed response? Gaussian LME on an asinh scale, where
%                         the zero mass is no longer present. Secondary.
%
% Contacts nest in electrode shafts and shafts in patients; both levels carry
% random intercepts. Linear-model denominator degrees of freedom use
% Satterthwaite so that inference is not silently carried out at the contact
% level. Region contrast families are FDR-controlled (Benjamini-Hochberg).
% Anything conditioned on the legacy selection rule is written to a separately
% named table and labelled secondary.

arguments
    T table
    cfg (1,1) struct
end

%% ------------------------------------------------------------------ Assemble
region = phg.cleanRegionLabels(T.NMM);
[regionNames, ~, regionIndex] = unique(region);
regionCounts = accumarray(regionIndex, 1);
[~, order] = sort(regionCounts, 'descend');
topCount = min(cfg.topRegionCount, numel(order));
topRegions = regionNames(order(1:topCount));
topRegions = setdiff(topRegions, "Other", 'stable');
regionGroup = region;
regionGroup(~ismember(regionGroup, topRegions)) = "Other";

ptID = categorical(string(T.PtID));
lead = phg.parseLeadLabel(T.Label);
leadID = categorical(string(T.PtID) + "/" + lead);

% "Other" is the pooled remainder of the sampled brain and is used as the
% reference level, so each named region is contrasted against the rest of the
% implant rather than against an arbitrary alphabetically-first parcel.
regionGroup = categorical(regionGroup);
otherFirst = ["Other", setdiff(string(categories(regionGroup))', "Other", ...
    'stable')];
regionGroup = reordercats(regionGroup, otherFirst);

signedRaw = T.RespAreaNet;
hasExcursion = signedRaw ~= 0;
isDilation = signedRaw > 0;
legacySelected = T.RespSig > cfg.legacySelectionThreshold;

isHippocampus = regionGroup == "Hippocampus";
contactClass = repmat("Extrahippocampal", height(T), 1);
contactClass(isHippocampus) = "Hippocampal";
contactClass = categorical(contactClass, ["Extrahippocampal", "Hippocampal"]);

% Robust standardisation for the magnitude model. asinh is used rather than log
% because the signed response is symmetric about zero and heavy-tailed in both
% directions. The scale is estimated from contacts that actually have an
% excursion, because the zero mass drives the pooled MAD to exactly zero.
excursionValues = signedRaw(hasExcursion);
center = median(excursionValues, 'omitnan');
scale = median(abs(excursionValues - center), 'omitnan');
if ~isfinite(scale) || scale <= eps
    scale = std(excursionValues, 'omitnan');
end
if ~isfinite(scale) || scale <= eps
    scale = 1;
end
signedCoupling = asinh((signedRaw - center) ./ scale);

modelTable = table(ptID, leadID, regionGroup, contactClass, hasExcursion, ...
    isDilation, signedRaw, signedCoupling, legacySelected, T.RespSig, ...
    'VariableNames', {'PtID', 'Lead', 'Region', 'ContactClass', ...
    'HasExcursion', 'IsDilation', 'SignedRaw', 'SignedCoupling', ...
    'LegacySelected', 'RespSig'});
modelTable = rmmissing(modelTable, 'DataVariables', 'SignedRaw');
excursionTable = modelTable(modelTable.HasExcursion, :);

nPatients = numel(categories(removecats(modelTable.PtID)));
nShafts = numel(categories(removecats(modelTable.Lead)));
fprintf(['[PHG] Model frame: %d contacts, %d patients, %d shafts; ', ...
    '%d contacts (%.1f%%) show a suprathreshold excursion.\n'], ...
    height(modelTable), nPatients, nShafts, height(excursionTable), ...
    100 * height(excursionTable) / height(modelTable));

%% ------------------------------------------------- Descriptive coverage
[group, patient, groupedRegion] = findgroups(string(modelTable.PtID), ...
    string(modelTable.Region));
nContacts = splitapply(@numel, modelTable.SignedRaw, group);
nExcursion = splitapply(@sum, modelTable.HasExcursion, group);
nLegacySelected = splitapply(@sum, modelTable.LegacySelected, group);
nDilation = splitapply(@(h, d) sum(h & d), modelTable.HasExcursion, ...
    modelTable.IsDilation, group);
medianSignedRaw = splitapply(@(x) median(x, 'omitnan'), ...
    modelTable.SignedRaw, group);
patientRegion = table(patient, groupedRegion, nContacts, nExcursion, ...
    nDilation, nLegacySelected, nExcursion ./ nContacts, medianSignedRaw, ...
    'VariableNames', {'patient_id', 'region', 'n_contacts', ...
    'n_excursion_contacts', 'n_dilation_contacts', 'n_legacy_selected', ...
    'excursion_proportion', 'median_signed_area'});
phg.writeTableAtomic(patientRegion, ...
    fullfile(cfg.tableDir, 'patient_region_summary.csv'));

[regionGroupId, summaryRegion] = findgroups(string(modelTable.Region));
coverage = splitapply(@numel, modelTable.SignedRaw, regionGroupId);
excursionCount = splitapply(@sum, modelTable.HasExcursion, regionGroupId);
dilationCount = splitapply(@(h, d) sum(h & d), modelTable.HasExcursion, ...
    modelTable.IsDilation, regionGroupId);
selectedCount = splitapply(@sum, modelTable.LegacySelected, regionGroupId);
patientCount = splitapply(@(x) numel(unique(x)), ...
    string(modelTable.PtID), regionGroupId);
shaftCount = splitapply(@(x) numel(unique(x)), ...
    string(modelTable.Lead), regionGroupId);
regionSummary = table(summaryRegion, coverage, excursionCount, ...
    dilationCount, excursionCount - dilationCount, selectedCount, ...
    patientCount, shaftCount, excursionCount ./ coverage, ...
    dilationCount ./ max(excursionCount, 1), 'VariableNames', {'region', ...
    'n_contacts', 'n_excursion_contacts', 'n_dilation_contacts', ...
    'n_constriction_contacts', 'n_legacy_selected', 'n_patients', ...
    'n_shafts', 'excursion_proportion', 'dilation_proportion'});
regionSummary = sortrows(regionSummary, 'n_contacts', 'descend');
phg.writeTableAtomic(regionSummary, ...
    fullfile(cfg.tableDir, 'region_coverage_summary.csv'));

%% ------------------------------ Part 1: prevalence of any excursion
prevalenceModel = fitglme(modelTable, ...
    'HasExcursion ~ Region + (1|PtID) + (1|PtID:Lead)', ...
    'Distribution', 'Binomial', 'Link', 'logit', 'FitMethod', 'Laplace');
prevalenceTable = localOddsTable(prevalenceModel, "Region_");
phg.writeTableAtomic(prevalenceTable, ...
    fullfile(cfg.tableDir, 'excursion_prevalence_glme.csv'));

%% ---------------- Part 2 (primary): direction given an excursion
% Hippocampal versus extrahippocampal is the prespecified single-degree-of-
% freedom contrast. It is estimated on every contact with an excursion, so it
% does not inherit the arbitrary 0.10 contiguity cut.
directionModel = fitglme(excursionTable, ...
    'IsDilation ~ ContactClass + (1|PtID) + (1|PtID:Lead)', ...
    'Distribution', 'Binomial', 'Link', 'logit', 'FitMethod', 'Laplace');
directionTable = localOddsTable(directionModel, "ContactClass");
phg.writeTableAtomic(directionTable, ...
    fullfile(cfg.tableDir, 'polarity_primary_glme.csv'));

% Region-resolved polarity is exploratory: several regions are completely
% separated (every excursion contact has the same sign), so their coefficients
% are not estimable and are flagged rather than reported as effects.
regionDirectionModel = fitglme(excursionTable, ...
    'IsDilation ~ Region + (1|PtID)', 'Distribution', 'Binomial', ...
    'Link', 'logit', 'FitMethod', 'Laplace');
regionDirectionTable = localOddsTable(regionDirectionModel, "Region_");
separated = false(height(regionDirectionTable), 1);
for k = 1:height(regionDirectionTable)
    name = regionDirectionTable.term(k);
    if ~startsWith(name, "Region_")
        continue
    end
    thisRegion = extractAfter(name, "Region_");
    values = excursionTable.IsDilation(string(excursionTable.Region) == thisRegion);
    separated(k) = numel(unique(values)) < 2;
end
regionDirectionTable.complete_separation = separated;
phg.writeTableAtomic(regionDirectionTable, ...
    fullfile(cfg.tableDir, 'polarity_by_region_exploratory.csv'));

%% ------------- Part 2 robustness: patient as the unit of replication
% A within-patient paired comparison removes every between-patient confound and
% uses the patient count, not the contact count, as the sample size.
patients = categories(removecats(excursionTable.PtID));
pairedHipp = nan(numel(patients), 1);
pairedExtra = nan(numel(patients), 1);
pairedHippN = nan(numel(patients), 1);
pairedExtraN = nan(numel(patients), 1);
for k = 1:numel(patients)
    rows = excursionTable.PtID == patients{k};
    hippRows = rows & excursionTable.ContactClass == "Hippocampal";
    extraRows = rows & excursionTable.ContactClass == "Extrahippocampal";
    pairedHippN(k) = sum(hippRows);
    pairedExtraN(k) = sum(extraRows);
    if any(hippRows)
        pairedHipp(k) = mean(excursionTable.IsDilation(hippRows));
    end
    if any(extraRows)
        pairedExtra(k) = mean(excursionTable.IsDilation(extraRows));
    end
end
complete = ~isnan(pairedHipp) & ~isnan(pairedExtra);
[pairedP, ~, pairedStats] = signrank(pairedHipp(complete), ...
    pairedExtra(complete));
pairedSigned = pairedStats.signedrank;

% Exact sign test on the same pairs. It discards effect magnitude, so a single
% discordant patient with a large reversal cannot dominate the statistic.
pairedDifference = pairedHipp(complete) - pairedExtra(complete);
nLower = sum(pairedDifference < 0);
nHigher = sum(pairedDifference > 0);
nTied = sum(pairedDifference == 0);
if nLower + nHigher > 0
    signTestP = 2 * min(binocdf(min(nLower, nHigher), nLower + nHigher, 0.5), 0.5);
    signTestP = min(signTestP, 1);
else
    signTestP = NaN;
end

pairedTable = table(string(patients), pairedHippN, pairedHipp, ...
    pairedExtraN, pairedExtra, pairedHipp - pairedExtra, ...
    'VariableNames', {'patient_id', 'n_hippocampal_excursion', ...
    'hippocampal_dilation_rate', 'n_extrahippocampal_excursion', ...
    'extrahippocampal_dilation_rate', 'dilation_rate_difference'});
phg.writeTableAtomic(pairedTable, ...
    fullfile(cfg.tableDir, 'polarity_within_patient_paired.csv'));
fprintf(['[PHG] Within-patient polarity: %d of %d informative patients ', ...
    'lower in hippocampus (%d tied); signed rank=%g, p=%.3g; ', ...
    'exact sign test p=%.3g\n'], nLower, nLower + nHigher, nTied, ...
    pairedSigned, pairedP, signTestP);

%% ------------------- Part 3 (secondary): magnitude given an excursion
magnitudeModel = fitlme(excursionTable, ...
    'SignedCoupling ~ Region + (1|PtID) + (1|PtID:Lead)', 'FitMethod', 'REML');
% fitlme reports residual denominator degrees of freedom by default, which
% would put the contact count in the denominator. Satterthwaite must be
% requested when the coefficients are extracted.
[~, ~, magnitudeStats] = fixedEffects(magnitudeModel, ...
    'DFMethod', 'satterthwaite');
termName = string(magnitudeStats.Name);
isRegionTerm = startsWith(termName, "Region_");
familyP = nan(height(magnitudeStats), 1);
familyP(isRegionTerm) = magnitudeStats.pValue(isRegionTerm);
magnitudeTable = table(termName, magnitudeStats.Estimate, magnitudeStats.SE, ...
    magnitudeStats.tStat, magnitudeStats.DF, magnitudeStats.pValue, ...
    phg.benjaminiHochberg(familyP), magnitudeStats.Lower, ...
    magnitudeStats.Upper, isRegionTerm, 'VariableNames', {'term', ...
    'estimate', 'standard_error', 't_statistic', 'satterthwaite_df', ...
    'p_value', 'fdr_q_value', 'ci95_low', 'ci95_high', 'in_fdr_family'});
phg.writeTableAtomic(magnitudeTable, ...
    fullfile(cfg.tableDir, 'signed_magnitude_mixed_model.csv'));
phg.writeText(fullfile(cfg.tableDir, 'signed_magnitude_mixed_model.txt'), ...
    evalc('disp(magnitudeModel)'), true);

%% ------------------- Marginal means on an interpretable scale
beta = magnitudeStats.Estimate;
V = magnitudeModel.CoefficientCovariance;
regionLevels = string(categories(removecats(excursionTable.Region)));
nLevels = numel(regionLevels);
emmEstimate = nan(nLevels, 1);
emmLow = nan(nLevels, 1);
emmHigh = nan(nLevels, 1);
emmDF = nan(nLevels, 1);
for k = 1:nLevels
    L = zeros(1, numel(beta));
    L(1) = 1;
    if regionLevels(k) == "Other"
        termDF = magnitudeStats.DF(1);
    else
        idx = find(termName == "Region_" + regionLevels(k), 1);
        if isempty(idx)
            continue
        end
        L(idx) = 1;
        termDF = magnitudeStats.DF(idx);
    end
    emmEstimate(k) = L * beta;
    standardError = sqrt(L * V * L');
    emmDF(k) = termDF;
    tCritical = tinv(0.975, termDF);
    emmLow(k) = emmEstimate(k) - tCritical * standardError;
    emmHigh(k) = emmEstimate(k) + tCritical * standardError;
end
backTransform = @(x) center + scale .* sinh(x);
marginalMeans = table(regionLevels, emmEstimate, emmLow, emmHigh, emmDF, ...
    backTransform(emmEstimate), backTransform(emmLow), ...
    backTransform(emmHigh), 'VariableNames', {'region', 'emm_asinh', ...
    'emm_asinh_ci95_low', 'emm_asinh_ci95_high', 'satterthwaite_df', ...
    'emm_signed_area', 'emm_signed_area_ci95_low', ...
    'emm_signed_area_ci95_high'});
marginalMeans = sortrows(marginalMeans, 'emm_asinh', 'descend');
phg.writeTableAtomic(marginalMeans, ...
    fullfile(cfg.tableDir, 'region_marginal_means.csv'));

%% ------------------------------------------------ Prespecified contrasts
contrastSpec = [
    "MTG middle temporal gyrus", "Hippocampus"
    "MTG middle temporal gyrus", "Other"
    "Hippocampus", "Other"
    "PrG precentral gyrus", "Other"
    ];
nContrast = size(contrastSpec, 1);
contrastName = strings(nContrast, 1);
contrastEstimate = nan(nContrast, 1);
contrastSE = nan(nContrast, 1);
contrastLow = nan(nContrast, 1);
contrastHigh = nan(nContrast, 1);
contrastP = nan(nContrast, 1);
contrastDF = nan(nContrast, 1);
for k = 1:nContrast
    L = zeros(1, numel(beta));
    positive = contrastSpec(k, 1);
    negative = contrastSpec(k, 2);
    if positive ~= "Other"
        L(termName == "Region_" + positive) = 1;
    end
    if negative ~= "Other"
        L(termName == "Region_" + negative) = -1;
    end
    contrastName(k) = positive + " minus " + negative;
    contrastEstimate(k) = L * beta;
    contrastSE(k) = sqrt(L * V * L');
    [pValue, ~, ~, denominatorDF] = coefTest(magnitudeModel, L, 0, ...
        'DFMethod', 'satterthwaite');
    contrastP(k) = pValue;
    contrastDF(k) = denominatorDF;
    tCrit = tinv(0.975, denominatorDF);
    contrastLow(k) = contrastEstimate(k) - tCrit * contrastSE(k);
    contrastHigh(k) = contrastEstimate(k) + tCrit * contrastSE(k);
end
contrastTable = table(contrastName, contrastEstimate, contrastSE, ...
    contrastLow, contrastHigh, contrastDF, contrastP, ...
    phg.benjaminiHochberg(contrastP), 'VariableNames', {'contrast', ...
    'estimate_asinh', 'standard_error', 'ci95_low', 'ci95_high', ...
    'satterthwaite_df', 'p_value', 'fdr_q_value'});
phg.writeTableAtomic(contrastTable, ...
    fullfile(cfg.tableDir, 'signed_magnitude_contrasts.csv'));

%% ------------------------------- Variance components and design summary
randomStats = covarianceParameters(magnitudeModel);
varianceRows = strings(0, 1);
varianceValue = [];
for k = 1:numel(randomStats)
    entry = randomStats{k};
    if istable(entry)
        for r = 1:height(entry)
            varianceRows(end+1, 1) = "magnitude_variance_" + ...
                string(entry.Group(r)); %#ok<AGROW>
            varianceValue(end+1, 1) = entry.Estimate(r); %#ok<AGROW>
        end
    end
end
designSummary = table([varianceRows; "n_contacts"; "n_excursion_contacts"; ...
    "n_patients"; "n_shafts"; "asinh_center"; "asinh_scale"; ...
    "paired_patients"; "paired_signed_rank"; "paired_p_value"; ...
    "paired_patients_lower_in_hippocampus"; ...
    "paired_patients_informative"; "paired_sign_test_p_value"], ...
    [varianceValue; height(modelTable); height(excursionTable); nPatients; ...
    nShafts; center; scale; sum(complete); pairedSigned; pairedP; ...
    nLower; nLower + nHigher; signTestP], ...
    'VariableNames', {'parameter', 'value'});
phg.writeTableAtomic(designSummary, ...
    fullfile(cfg.tableDir, 'analysis_design_summary.csv'));

%% --------------------- Sensitivity: does polarity survive the threshold?
% Threshold 0 is the hurdle model itself. The remaining rows re-estimate the
% dissociation under progressively stricter versions of the legacy contiguity
% cut. A conclusion that moves with the cut would be an artefact of the cut.
thresholds = [0; 0.05; 0.10; 0.15; 0.20; 0.30];
sensitivityN = nan(size(thresholds));
sensitivityHippN = nan(size(thresholds));
sensitivityOdds = nan(size(thresholds));
sensitivityLow = nan(size(thresholds));
sensitivityHigh = nan(size(thresholds));
sensitivityP = nan(size(thresholds));
hippDilationRate = nan(size(thresholds));
extraDilationRate = nan(size(thresholds));
for k = 1:numel(thresholds)
    frame = modelTable(modelTable.HasExcursion & ...
        modelTable.RespSig > thresholds(k), :);
    sensitivityN(k) = height(frame);
    isHipp = frame.ContactClass == "Hippocampal";
    sensitivityHippN(k) = sum(isHipp);
    hippDilationRate(k) = mean(frame.IsDilation(isHipp));
    extraDilationRate(k) = mean(frame.IsDilation(~isHipp));
    if sensitivityHippN(k) < 5 || numel(unique(frame.IsDilation)) < 2
        continue
    end
    try
        model = fitglme(frame, 'IsDilation ~ ContactClass + (1|PtID)', ...
            'Distribution', 'Binomial', 'Link', 'logit', ...
            'FitMethod', 'Laplace');
        idx = find(startsWith(string(model.Coefficients.Name), ...
            "ContactClass"), 1);
        if isempty(idx)
            continue
        end
        sensitivityOdds(k) = exp(model.Coefficients.Estimate(idx));
        sensitivityLow(k) = exp(model.Coefficients.Lower(idx));
        sensitivityHigh(k) = exp(model.Coefficients.Upper(idx));
        sensitivityP(k) = model.Coefficients.pValue(idx);
    catch modelError
        warning('phg:PolaritySensitivity', ...
            'Threshold %.2f did not converge: %s', thresholds(k), ...
            modelError.message);
    end
end
sensitivityTable = table(thresholds, sensitivityN, sensitivityHippN, ...
    hippDilationRate, extraDilationRate, sensitivityOdds, sensitivityLow, ...
    sensitivityHigh, sensitivityP, 'VariableNames', {'selection_threshold', ...
    'n_contacts', 'n_hippocampal_contacts', 'hippocampal_dilation_rate', ...
    'extrahippocampal_dilation_rate', 'odds_ratio_hippocampal', ...
    'odds_ratio_ci95_low', 'odds_ratio_ci95_high', 'p_value'});
phg.writeTableAtomic(sensitivityTable, ...
    fullfile(cfg.tableDir, 'polarity_threshold_sensitivity.csv'));

outputs = struct('patientRegion', patientRegion, ...
    'regionSummary', regionSummary, 'prevalenceModel', prevalenceModel, ...
    'prevalenceTable', prevalenceTable, 'directionModel', directionModel, ...
    'directionTable', directionTable, ...
    'regionDirectionTable', regionDirectionTable, ...
    'pairedTable', pairedTable, 'pairedP', pairedP, ...
    'magnitudeModel', magnitudeModel, 'magnitudeTable', magnitudeTable, ...
    'marginalMeans', marginalMeans, 'contrasts', contrastTable, ...
    'designSummary', designSummary, 'polaritySensitivity', sensitivityTable);

fprintf(['[PHG] Hurdle analysis complete: prevalence, direction ', ...
    '(primary), and magnitude models written.\n']);
end

% -------------------------------------------------------------------------
function outTable = localOddsTable(model, familyPrefix)
%LOCALODDSTABLE Coefficient table on the odds-ratio scale with FDR control.

stats = model.Coefficients;
name = string(stats.Name);
inFamily = startsWith(name, familyPrefix);
familyP = nan(numel(name), 1);
familyP(inFamily) = stats.pValue(inFamily);
outTable = table(name, stats.Estimate, stats.SE, exp(stats.Estimate), ...
    exp(stats.Lower), exp(stats.Upper), stats.DF, stats.pValue, ...
    phg.benjaminiHochberg(familyP), inFamily, 'VariableNames', {'term', ...
    'log_odds', 'standard_error', 'odds_ratio', 'odds_ratio_ci95_low', ...
    'odds_ratio_ci95_high', 'df', 'p_value', 'fdr_q_value', 'in_fdr_family'});
end
