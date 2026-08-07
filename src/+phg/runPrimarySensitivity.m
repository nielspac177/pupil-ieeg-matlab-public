function outputs = runPrimarySensitivity(channel, cfg)
%RUNPRIMARYSENSITIVITY Does the polarity effect survive adjustment?
%
%   The primary model contains one fixed effect. A reviewer will ask whether
%   the effect is really about the hippocampus or about something correlated
%   with where hippocampal electrodes go and how they are recorded. Three
%   adjustments are available without new data, and each is fitted as a
%   separate model rather than all at once, because with 285 excursion
%   contacts a four-covariate model would be fitting noise.
%
%     condition   pupil behaves differently at rest and in a task, and the
%                 cohort is not balanced across the three conditions
%     hemisphere  mesial temporal implants are usually lateralised toward the
%                 suspected focus, so side and pathology travel together
%     peak count  a contact with 313 detected high-gamma peaks has a noisier
%                 peri-peak estimate than one with 1000; log peak count enters
%                 as a proxy for how well each contact's response is measured
%
%   A fourth question is asked instead of adjusted for: among hippocampal
%   contacts alone, does the constriction preference differ between left and
%   right? A large asymmetry would point at pathology rather than physiology.
%
%   Every model here is secondary. The primary estimate is the unadjusted one,
%   declared in advance; these say whether it moves.

arguments
    channel table
    cfg (1,1) struct
end

region = phg.cleanRegionLabels(channel.NMM);
label = string(channel.Label);

frame = table;
frame.PtID = categorical(string(channel.PtID));
frame.Lead = categorical(phg.parseLeadLabel(label));
frame.IsDilation = double(channel.RespAreaNet > 0);
frame.HasExcursion = channel.RespAreaNet ~= 0;
frame.ContactClass = categorical( ...
    repmat("Extrahippocampal", height(channel), 1));
frame.ContactClass(region == "Hippocampus") = "Hippocampal";
frame.ContactClass = reordercats(frame.ContactClass, ...
    {'Extrahippocampal', 'Hippocampal'});

% Hemisphere from the contact label prefix. It agrees with the sign of the MNI
% x coordinate on every contact that carries a prefix, so the label is used and
% the coordinate is the check rather than the source.
side = repmat("Unknown", height(channel), 1);
side(startsWith(label, "L")) = "Left";
side(startsWith(label, "R")) = "Right";
coordinateSide = repmat("Left", height(channel), 1);
coordinateSide(channel.XYZMNI(:,1) > 0) = "Right";
known = side ~= "Unknown";
agreement = mean(side(known) == coordinateSide(known));
if agreement < 0.95
    warning('phg:lateralityDisagreement', ...
        ['Label-derived hemisphere agrees with the MNI x coordinate on only ' ...
         '%.0f%% of contacts; not using it.'], 100*agreement);
    side(:) = "Unknown";
end
frame.Hemisphere = categorical(side);

taskMap = fullfile(cfg.repoRoot, 'config', 'patient_task_map.tsv');
if isfile(taskMap)
    map = readtable(taskMap, 'FileType', 'text', 'Delimiter', '\t', ...
        'TextType', 'string');
    condition = repmat("Unknown", height(channel), 1);
    [found, where] = ismember(string(channel.PtID), map.patient_id);
    condition(found) = map.task(where(found));
    frame.Condition = categorical(condition);
end

% The number of detected peaks is the sample size behind each contact's
% peri-peak average; it is skewed, so it enters on a log scale and centred.
peakCount = double(channel.N);
frame.LogPeaks = log(peakCount) - mean(log(peakCount), 'omitnan');

excursion = frame(frame.HasExcursion, :);
excursion.PtID = removecats(excursion.PtID);

models = { ...
    'unadjusted',        'IsDilation ~ ContactClass + (1|PtID) + (1|PtID:Lead)'; ...
    'plus condition',    'IsDilation ~ ContactClass + Condition + (1|PtID) + (1|PtID:Lead)'; ...
    'plus hemisphere',   'IsDilation ~ ContactClass + Hemisphere + (1|PtID) + (1|PtID:Lead)'; ...
    'plus log peaks',    'IsDilation ~ ContactClass + LogPeaks + (1|PtID) + (1|PtID:Lead)'};

rows = cell(0, 1);
for k = 1:size(models, 1)
    formula = models{k, 2};
    if contains(formula, 'Condition') && ~ismember('Condition', ...
            excursion.Properties.VariableNames)
        continue
    end
    try
        model = fitglme(excursion, formula, 'Distribution', 'Binomial', ...
            'Link', 'logit', 'FitMethod', 'Laplace');
    catch modelError
        warning('phg:sensitivityFitFailed', '%s: %s', models{k,1}, ...
            modelError.message);
        continue
    end
    coefficients = model.Coefficients;
    index = find(startsWith(string(coefficients.Name), "ContactClass"), 1);
    if isempty(index)
        continue
    end
    rows{end+1, 1} = table(string(models{k,1}), height(excursion), ...
        exp(coefficients.Estimate(index)), exp(coefficients.Lower(index)), ...
        exp(coefficients.Upper(index)), coefficients.pValue(index), ...
        'VariableNames', {'model', 'n_contacts', 'odds_ratio', ...
        'odds_ratio_ci95_low', 'odds_ratio_ci95_high', 'p_value'}); %#ok<AGROW>
end
resultTable = vertcat(rows{:});
phg.writeTableAtomic(resultTable, ...
    fullfile(cfg.tableDir, 'polarity_adjusted_models.csv'));

% Laterality within the hippocampus: asked, not adjusted for.
hippocampal = excursion(excursion.ContactClass == "Hippocampal" & ...
    excursion.Hemisphere ~= "Unknown", :);
hippocampal.Hemisphere = removecats(hippocampal.Hemisphere);
hippocampal.PtID = removecats(hippocampal.PtID);
lateralityTable = table.empty;
if numel(categories(hippocampal.Hemisphere)) == 2 && ...
        numel(unique(hippocampal.IsDilation)) == 2
    lateralityModel = fitglme(hippocampal, ...
        'IsDilation ~ Hemisphere + (1|PtID)', 'Distribution', 'Binomial', ...
        'Link', 'logit', 'FitMethod', 'Laplace');
    coefficients = lateralityModel.Coefficients;
    index = find(startsWith(string(coefficients.Name), "Hemisphere"), 1);
    lateralityTable = table(string(coefficients.Name(index)), ...
        height(hippocampal), exp(coefficients.Estimate(index)), ...
        exp(coefficients.Lower(index)), exp(coefficients.Upper(index)), ...
        coefficients.pValue(index), ...
        'VariableNames', {'term', 'n_hippocampal_contacts', 'odds_ratio', ...
        'odds_ratio_ci95_low', 'odds_ratio_ci95_high', 'p_value'});
else
    lateralityTable = table("not estimable", height(hippocampal), NaN, ...
        NaN, NaN, NaN, 'VariableNames', {'term', 'n_hippocampal_contacts', ...
        'odds_ratio', 'odds_ratio_ci95_low', 'odds_ratio_ci95_high', ...
        'p_value'});
end
phg.writeTableAtomic(lateralityTable, ...
    fullfile(cfg.tableDir, 'hippocampal_laterality.csv'));

outputs = struct('adjustedTable', resultTable, ...
    'lateralityTable', lateralityTable);

for k = 1:height(resultTable)
    fprintf('[PHG] Polarity, %-16s OR %.3f [%.3f, %.3f], P = %.3g\n', ...
        resultTable.model(k), resultTable.odds_ratio(k), ...
        resultTable.odds_ratio_ci95_low(k), ...
        resultTable.odds_ratio_ci95_high(k), resultTable.p_value(k));
end
if ~isnan(lateralityTable.odds_ratio(1))
    fprintf(['[PHG] Hippocampal laterality: %s OR %.2f [%.2f, %.2f], ' ...
        'P = %.3g on %d contacts.\n'], lateralityTable.term(1), ...
        lateralityTable.odds_ratio(1), lateralityTable.odds_ratio_ci95_low(1), ...
        lateralityTable.odds_ratio_ci95_high(1), lateralityTable.p_value(1), ...
        lateralityTable.n_hippocampal_contacts(1));
end
end
