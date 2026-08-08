function outputs = runReplicationAnalyses(derived, cfg, options)
%RUNREPLICATIONANALYSES H1, H3 and the prespecified decision rule.
%
% The decision rule is applied mechanically, from the numbers, exactly as it was
% written in docs/replication_plan_ebrains.md section 2 before the data arrived:
%
%   replication              OR < 1 and the 95% interval excludes 1
%   directionally consistent OR < 1, interval includes 1, and the point
%     and underpowered       estimate lies inside the discovery interval
%                            [0.021, 0.207]
%   inconclusive             OR < 1, interval includes 1, point estimate
%                            outside the discovery interval
%   failure                  OR >= 1, or an interval excluding 1 in the
%                            opposite direction
%
% The verdict is computed here rather than written by hand so that it cannot
% drift when the numbers change. It is exported to a CSV like every other
% quantity in this project.

arguments
    derived table
    cfg (1,1) struct
    options.label (1,1) string = "primary"
    options.writeTables (1,1) logical = true
end

frame = localModelFrame(derived);
excursion = frame(frame.HasExcursion, :);

nHippocampalExcursion = sum(excursion.ContactClass == "Hippocampal");
fprintf(['[PHG] Replication frame (%s): %d contacts, %d subjects, ', ...
    '%d shafts; %d excursion contacts (%d hippocampal).\n'], ...
    options.label, height(frame), numel(categories(removecats(frame.PtID))), ...
    numel(categories(removecats(frame.Lead))), height(excursion), ...
    nHippocampalExcursion);

outputs = struct('frame', frame, 'excursion', excursion);

% Amendment section 7.2: below this the primary model is not estimable, and
% that is reported as such rather than as a null.
if nHippocampalExcursion < cfg.replication.minHippocampalExcursion || ...
        numel(unique(excursion.ContactClass)) < 2
    verdict = localVerdictTable(options.label, NaN, NaN, NaN, NaN, ...
        height(excursion), nHippocampalExcursion, NaN, NaN, ...
        "not_estimable", cfg);
    if options.writeTables
        phg.writeTableAtomic(verdict, fullfile(cfg.tableDir, ...
            "replication_h1_decision_" + options.label + ".csv"));
    end
    % Callers read these unconditionally; a not-estimable result must still
    % produce a well-formed output rather than an error three steps later.
    outputs.verdict = verdict;
    outputs.coefficients = table();
    outputs.prevalence = table();
    outputs.adjusted = table();
    outputs.coverage = table();
    outputs.overlapGap = NaN;
    warning('phg:ReplicationNotEstimable', ...
        ['Only %d hippocampal excursion contacts; H1 is reported as not ', ...
        'estimable.'], nHippocampalExcursion);
    return
end

%% ------------------------------------------- H1 primary, confirmatory
model = fitglme(excursion, ...
    'IsDilation ~ ContactClass + (1|PtID) + (1|PtID:Lead)', ...
    'Distribution', 'Binomial', 'Link', 'logit', 'FitMethod', 'Laplace');
coefficients = localOddsTable(model, "ContactClass");
row = coefficients(startsWith(coefficients.term, "ContactClass"), :);

oddsRatio = row.odds_ratio(1);
ciLow = row.odds_ratio_ci95_low(1);
ciHigh = row.odds_ratio_ci95_high(1);
pValue = row.p_value(1);

hippocampalRate = mean(excursion.IsDilation( ...
    excursion.ContactClass == "Hippocampal"));
extrahippocampalRate = mean(excursion.IsDilation( ...
    excursion.ContactClass == "Extrahippocampal"));

decision = phg.replicationVerdict(oddsRatio, ciLow, ciHigh, ...
    cfg.replication.discoveryInterval);
verdict = localVerdictTable(options.label, oddsRatio, ciLow, ciHigh, ...
    pValue, height(excursion), nHippocampalExcursion, hippocampalRate, ...
    extrahippocampalRate, decision, cfg);

fprintf(['[PHG] H1 (%s): OR = %.3f [%.3f, %.3f], P = %.3g; ', ...
    '%.0f%% vs %.0f%% dilation-linked -> %s\n'], options.label, oddsRatio, ...
    ciLow, ciHigh, pValue, 100 * hippocampalRate, ...
    100 * extrahippocampalRate, decision);

%% ---------------------------- H1 adjusted for peri-peak event overlap
% Prespecified trigger: adjust only if the classes differ by more than five
% percentage points in how much of the peri-peak window sits near a task event.
overlapGap = abs(mean(excursion.EventOverlap( ...
    excursion.ContactClass == "Hippocampal"), 'omitnan') - ...
    mean(excursion.EventOverlap( ...
    excursion.ContactClass == "Extrahippocampal"), 'omitnan'));
adjusted = table();
if isfinite(overlapGap) && overlapGap > cfg.replication.overlapAdjustmentGap
    adjustedModel = fitglme(excursion, ...
        ['IsDilation ~ ContactClass + EventOverlap + (1|PtID) + ' ...
        '(1|PtID:Lead)'], 'Distribution', 'Binomial', 'Link', 'logit', ...
        'FitMethod', 'Laplace');
    adjusted = localOddsTable(adjustedModel, "ContactClass");
end

%% ----------------------------------------------- H3, exploratory
prevalence = table();
try
    prevalenceModel = fitglme(frame, ...
        'HasExcursion ~ ContactClass + (1|PtID) + (1|PtID:Lead)', ...
        'Distribution', 'Binomial', 'Link', 'logit', 'FitMethod', 'Laplace');
    prevalence = localOddsTable(prevalenceModel, "ContactClass");
catch prevalenceError
    warning('phg:ReplicationPrevalence', ...
        'H3 prevalence model did not converge: %s', prevalenceError.message);
end

%% ------------------------------------------- Descriptive coverage
[groupId, subject, region] = findgroups(string(frame.PtID), ...
    string(frame.Region));
coverage = table(subject, region, splitapply(@numel, frame.SignedRaw, groupId), ...
    splitapply(@sum, frame.HasExcursion, groupId), ...
    splitapply(@(h, d) sum(h & d), frame.HasExcursion, frame.IsDilation, ...
    groupId), 'VariableNames', {'subject_id', 'region', 'n_contacts', ...
    'n_excursion_contacts', 'n_dilation_contacts'});

if options.writeTables
    suffix = "_" + options.label;
    phg.writeTableAtomic(coefficients, fullfile(cfg.tableDir, ...
        "replication_h1_glme" + suffix + ".csv"));
    phg.writeTableAtomic(verdict, fullfile(cfg.tableDir, ...
        "replication_h1_decision" + suffix + ".csv"));
    phg.writeTableAtomic(coverage, fullfile(cfg.tableDir, ...
        "replication_coverage" + suffix + ".csv"));
    if ~isempty(prevalence)
        phg.writeTableAtomic(prevalence, fullfile(cfg.tableDir, ...
            "replication_h3_prevalence" + suffix + ".csv"));
    end
    if ~isempty(adjusted)
        phg.writeTableAtomic(adjusted, fullfile(cfg.tableDir, ...
            "replication_h1_overlap_adjusted" + suffix + ".csv"));
    end
end

outputs.model = model;
outputs.coefficients = coefficients;
outputs.verdict = verdict;
outputs.prevalence = prevalence;
outputs.adjusted = adjusted;
outputs.coverage = coverage;
outputs.overlapGap = overlapGap;
end

% -------------------------------------------------------------------------
function frame = localModelFrame(derived)
%LOCALMODELFRAME Assemble the modelling frame in the discovery's own terms.

region = phg.cleanRegionLabels(derived.NMM);
ptID = categorical(string(derived.PtID));
leadID = categorical(string(derived.PtID) + "/" + string(derived.Shaft));

signedRaw = derived.RespAreaNet;
hasExcursion = signedRaw ~= 0;
isDilation = signedRaw > 0;

% Hippocampus proper. Parahippocampal gyrus and entorhinal cortex are distinct
% parcels and belong to the reference class, matching loadDerivedTables.m.
isHippocampus = lower(region) == "hippocampus";
contactClass = repmat("Extrahippocampal", height(derived), 1);
contactClass(isHippocampus) = "Hippocampal";
contactClass = categorical(contactClass, ...
    ["Extrahippocampal", "Hippocampal"]);

frame = table(ptID, leadID, categorical(region), contactClass, hasExcursion, ...
    isDilation, signedRaw, derived.RespSig, derived.EventOverlap, ...
    'VariableNames', {'PtID', 'Lead', 'Region', 'ContactClass', ...
    'HasExcursion', 'IsDilation', 'SignedRaw', 'RespSig', 'EventOverlap'});
frame = rmmissing(frame, 'DataVariables', 'SignedRaw');
end

% -------------------------------------------------------------------------
function verdict = localVerdictTable(label, oddsRatio, ciLow, ciHigh, ...
    pValue, nExcursion, nHippocampal, hippocampalRate, extrahippocampalRate, ...
    decision, cfg)
%LOCALVERDICTTABLE One row carrying the result and the rule that judged it.

discovery = cfg.replication.discoveryInterval;
verdict = table(string(label), oddsRatio, ciLow, ciHigh, pValue, ...
    nExcursion, nHippocampal, hippocampalRate, extrahippocampalRate, ...
    discovery(1), discovery(2), string(decision), ...
    'VariableNames', {'analysis', 'odds_ratio_hippocampal', ...
    'odds_ratio_ci95_low', 'odds_ratio_ci95_high', 'p_value', ...
    'n_excursion_contacts', 'n_hippocampal_excursion_contacts', ...
    'hippocampal_dilation_rate', 'extrahippocampal_dilation_rate', ...
    'discovery_ci_low', 'discovery_ci_high', 'prespecified_verdict'});
end

% -------------------------------------------------------------------------
function outTable = localOddsTable(model, familyPrefix)
%LOCALODDSTABLE Coefficients on the odds-ratio scale.

stats = model.Coefficients;
name = string(stats.Name);
inFamily = startsWith(name, familyPrefix);
outTable = table(name, stats.Estimate, stats.SE, exp(stats.Estimate), ...
    exp(stats.Lower), exp(stats.Upper), stats.DF, stats.pValue, inFamily, ...
    'VariableNames', {'term', 'log_odds', 'standard_error', 'odds_ratio', ...
    'odds_ratio_ci95_low', 'odds_ratio_ci95_high', 'df', 'p_value', ...
    'in_family'});
end
