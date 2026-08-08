function results = assembleReplicationResults(cfg, analyses, ripple)
%ASSEMBLEREPLICATIONRESULTS Collect the replication into one reportable object.
%   Everything the manuscript quotes comes from the tables written here, so no
%   replication number is ever typed into the text by hand.

arguments
    cfg (1,1) struct
    analyses struct
    ripple (1,1) struct
end

labels = string({analyses.label});
primaryIndex = find(labels == "primary", 1);
if isempty(primaryIndex)
    primaryIndex = 1;
end
primary = analyses(primaryIndex);

%% ------------------------------------------------------------- forest
prettyName = containers.Map( ...
    {'primary', 'raw_pupil', 'rule_b', 'rule_c'}, ...
    {'Primary: gaze-regressed, rule A', 'Sensitivity: raw pupil', ...
     'Sensitivity: rule B (±5 s clear)', 'Sensitivity: rule C (−20/+5 s clear)'});

label = strings(numel(analyses), 1);
oddsRatio = nan(numel(analyses), 1);
ciLow = nan(numel(analyses), 1);
ciHigh = nan(numel(analyses), 1);
pValue = nan(numel(analyses), 1);
nExcursion = nan(numel(analyses), 1);
nHippocampal = nan(numel(analyses), 1);
verdicts = strings(numel(analyses), 1);
isPrimary = false(numel(analyses), 1);

for k = 1:numel(analyses)
    key = char(analyses(k).label);
    if prettyName.isKey(key)
        label(k) = string(prettyName(key));
    else
        label(k) = string(key);
    end
    verdict = analyses(k).verdict;
    oddsRatio(k) = verdict.odds_ratio_hippocampal(1);
    ciLow(k) = verdict.odds_ratio_ci95_low(1);
    ciHigh(k) = verdict.odds_ratio_ci95_high(1);
    pValue(k) = verdict.p_value(1);
    nExcursion(k) = verdict.n_excursion_contacts(1);
    nHippocampal(k) = verdict.n_hippocampal_excursion_contacts(1);
    verdicts(k) = verdict.prespecified_verdict(1);
    isPrimary(k) = k == primaryIndex;
end

% Separation flag. A confidence interval spanning many orders of magnitude is
% not a wide estimate, it is a model that did not converge on data where one
% class has no variation left -- the maximum-likelihood estimate runs to
% infinity and the optimiser stops wherever it stops. The discovery analysis
% already uses this convention, reporting a hippocampal laterality odds ratio
% of 6.81 [0.02, 2294] as "not estimable" rather than as a null.
%
% This does not touch the H1 decision rule, which applies to the primary
% analysis alone. It marks sensitivity analyses whose numbers cannot be read as
% effect sizes. The threshold was set after seeing that one sensitivity fit
% returned an interval nineteen orders of magnitude wide, and it is recorded
% here rather than applied silently.
intervalDecades = log10(ciHigh ./ ciLow);
separated = isfinite(intervalDecades) & intervalDecades > 4;
verdicts(separated) = "not_estimable_separation";

forest = table(label, oddsRatio, ciLow, ciHigh, pValue, nExcursion, ...
    nHippocampal, intervalDecades, separated, verdicts, isPrimary, ...
    'VariableNames', {'label', 'odds_ratio', 'ci_low', 'ci_high', ...
    'p_value', 'n_excursion_contacts', 'n_hippocampal_excursion_contacts', ...
    'ci_width_decades', 'separation_suspected', 'prespecified_verdict', ...
    'is_primary'});
phg.writeTableAtomic(forest, ...
    fullfile(cfg.tableDir, 'replication_h1_all_configurations.csv'));

%% ------------------------------------------------------ per participant
excursion = primary.excursion;
subjects = categories(removecats(excursion.PtID));
hippocampalRate = nan(numel(subjects), 1);
extrahippocampalRate = nan(numel(subjects), 1);
hippocampalCount = zeros(numel(subjects), 1);
extrahippocampalCount = zeros(numel(subjects), 1);
for k = 1:numel(subjects)
    rows = excursion.PtID == subjects{k};
    hippocampal = rows & excursion.ContactClass == "Hippocampal";
    extrahippocampal = rows & excursion.ContactClass == "Extrahippocampal";
    hippocampalCount(k) = sum(hippocampal);
    extrahippocampalCount(k) = sum(extrahippocampal);
    if any(hippocampal)
        hippocampalRate(k) = mean(excursion.IsDilation(hippocampal));
    end
    if any(extrahippocampal)
        extrahippocampalRate(k) = mean(excursion.IsDilation(extrahippocampal));
    end
end
perSubject = table(string(subjects), hippocampalCount, hippocampalRate, ...
    extrahippocampalCount, extrahippocampalRate, 'VariableNames', ...
    {'subject_id', 'n_hippocampal_excursion', 'hippocampal_rate', ...
    'n_extrahippocampal_excursion', 'extrahippocampal_rate'});
phg.writeTableAtomic(perSubject, ...
    fullfile(cfg.tableDir, 'replication_per_subject_rates.csv'));

%% ------------------------------------------------------- cohort rates
hippocampalRows = excursion.ContactClass == "Hippocampal";
cohortRates = struct;
cohortRates.labels = {'Extrahippocampal', 'Hippocampal'};
cohortRates.values = [mean(excursion.IsDilation(~hippocampalRows)); ...
    mean(excursion.IsDilation(hippocampalRows))];
cohortRates.counts = [sum(~hippocampalRows); sum(hippocampalRows)];

%% ------------------------------------------------------------- ripple
rippleSummary = struct('observed', NaN, 'null', [], 'pValue', NaN, ...
    'nRipples', 0);
if ~isempty(ripple.test) && height(ripple.test) > 0
    rippleSummary.observed = ripple.test.observed_constriction_fraction(1);
    rippleSummary.pValue = ripple.test.p_value_one_sided(1);
    rippleSummary.nRipples = ripple.test.n_ripples(1);
    % The pooled circular-shift surrogates themselves, not a per-run summary of
    % them: the figure has to show the distribution the p-value came from.
    if isfield(ripple, 'nullDistribution')
        rippleSummary.null = ripple.nullDistribution( ...
            isfinite(ripple.nullDistribution));
    end
end

%% --------------------------------------------------------- headline
headline = table(forest.odds_ratio(primaryIndex), ...
    forest.ci_low(primaryIndex), forest.ci_high(primaryIndex), ...
    forest.p_value(primaryIndex), forest.n_excursion_contacts(primaryIndex), ...
    forest.n_hippocampal_excursion_contacts(primaryIndex), ...
    cohortRates.values(2), cohortRates.values(1), ...
    numel(subjects), sum(isfinite(hippocampalRate)), ...
    forest.prespecified_verdict(primaryIndex), ...
    'VariableNames', {'odds_ratio', 'ci95_low', 'ci95_high', 'p_value', ...
    'n_excursion_contacts', 'n_hippocampal_excursion_contacts', ...
    'hippocampal_dilation_rate', 'extrahippocampal_dilation_rate', ...
    'n_subjects', 'n_subjects_with_hippocampal_coverage', ...
    'prespecified_verdict'});
phg.writeTableAtomic(headline, ...
    fullfile(cfg.tableDir, 'replication_headline.csv'));

results = struct('forest', forest, 'perSubject', perSubject, ...
    'cohortRates', cohortRates, 'ripple', rippleSummary, ...
    'headline', headline, 'verdict', forest.prespecified_verdict(primaryIndex));
end
