function outputs = runContinuousRobustness(channel, cfg)
%RUNCONTINUOUSROBUSTNESS Robustness of the primary, on the primary's own scale.
%
% The adjustments a reader asks for -- does this survive covariates, does it
% survive dropping the most influential region, does it survive holding out any
% one patient -- were all answered on the dichotomised secondary, and reported
% as odds ratios. That left the paper defending its primary result with the
% robustness of a different analysis, and describing electrophysiology in an
% effect measure built for comparing event rates between groups of people.
%
% The same questions are asked here of the model that is actually primary: a
% continuous signed coupling amplitude, every contact, no gate. Coefficients
% are reported in standard deviations of the outcome, which is the scale the
% rest of the paper uses and the one this literature reports.
%
% Three checks, matching the three the secondary already carried:
%
%   * Covariate adjustment. Recording condition, hemisphere, and the log
%     number of detected high-gamma peaks, the last standing in for how
%     precisely each contact's peri-peak average is estimated.
%   * Reference-class composition. Middle temporal gyrus supplies the most
%     contacts and is the most strongly dilation-coupled region, so it is the
%     obvious candidate for carrying the contrast on its own; the model is
%     refitted with it removed from the comparison class.
%   * Hemisphere restriction, since a contrast between one deep structure and
%     the rest of an implant could in principle be a laterality effect.

arguments
    channel table
    cfg (1,1) struct
end

warningState = warning('off', 'all');
cleanup = onCleanup(@() warning(warningState));

region = phg.cleanRegionLabels(channel.NMM);
label = string(channel.Label);

scale = median(abs(channel.FitHeight), 'omitnan');
if ~isfinite(scale) || scale <= eps
    scale = 1;
end

frame = table;
frame.PtID = categorical(string(channel.PtID));
frame.Lead = categorical(string(channel.PtID) + "/" + ...
    phg.parseLeadLabel(label));
frame.Coupling = asinh(channel.FitHeight ./ scale);
frame.ContactClass = categorical( ...
    repmat("Extrahippocampal", height(channel), 1));
frame.ContactClass(region == "Hippocampus") = "Hippocampal";
frame.ContactClass = reordercats(frame.ContactClass, ...
    {'Extrahippocampal', 'Hippocampal'});
frame.Region = region;

% Hemisphere from the label prefix, checked against the MNI x coordinate, as
% in runPrimarySensitivity. The label is the source and the coordinate the
% check, because a label can be missing but is not usually wrong.
side = repmat("Unknown", height(channel), 1);
side(startsWith(label, "L")) = "Left";
side(startsWith(label, "R")) = "Right";
coordinateSide = repmat("Left", height(channel), 1);
coordinateSide(channel.XYZMNI(:, 1) > 0) = "Right";
known = side ~= "Unknown";
if mean(side(known) == coordinateSide(known)) < 0.95
    side(:) = "Unknown";
end
frame.Hemisphere = categorical(side);

taskMap = fullfile(cfg.repoRoot, 'config', 'patient_task_map.tsv');
hasCondition = false;
if isfile(taskMap)
    map = readtable(taskMap, 'FileType', 'text', 'Delimiter', '\t', ...
        'TextType', 'string');
    condition = repmat("Unknown", height(channel), 1);
    [found, where] = ismember(string(channel.PtID), map.patient_id);
    condition(found) = map.task(where(found));
    frame.Condition = categorical(condition);
    hasCondition = numel(categories(removecats(frame.Condition))) > 1;
end

peakCount = double(channel.N);
frame.LogPeaks = log(peakCount) - mean(log(peakCount), 'omitnan');
frame = frame(isfinite(frame.Coupling), :);

specification = {
    "unadjusted", frame, ...
        'Coupling ~ ContactClass + (1|PtID) + (1|PtID:Lead)'
    "plus hemisphere", frame(frame.Hemisphere ~= "Unknown", :), ...
        'Coupling ~ ContactClass + Hemisphere + (1|PtID) + (1|PtID:Lead)'
    "plus log peaks", frame, ...
        'Coupling ~ ContactClass + LogPeaks + (1|PtID) + (1|PtID:Lead)'
    "excluding middle temporal gyrus", ...
        frame(~startsWith(frame.Region, "MTG") | ...
              frame.ContactClass == "Hippocampal", :), ...
        'Coupling ~ ContactClass + (1|PtID) + (1|PtID:Lead)'
    };
if hasCondition
    specification(end + 1, :) = {"plus condition", frame, ...
        'Coupling ~ ContactClass + Condition + (1|PtID) + (1|PtID:Lead)'};
end

name = strings(0, 1);
nContacts = []; estimate = []; low = []; high = []; pValue = []; sdUnits = [];
for k = 1:size(specification, 1)
    subset = specification{k, 2};
    % Subsetting leaves categorical levels behind with no observations, and an
    % empty level still contributes a column to the design matrix, which is
    % then rank deficient. Drop unused levels on every categorical column.
    for column = string(subset.Properties.VariableNames)
        if iscategorical(subset.(column))
            subset.(column) = removecats(subset.(column));
        end
    end
    if height(subset) < 50 || numel(unique(subset.ContactClass)) < 2
        continue
    end
    model = fitlme(subset, specification{k, 3}, 'FitMethod', 'REML');
    [~, ~, stats] = fixedEffects(model, 'DFMethod', 'satterthwaite');
    index = find(startsWith(string(stats.Name), "ContactClass"), 1);
    if isempty(index)
        continue
    end
    name(end + 1, 1) = specification{k, 1}; %#ok<AGROW>
    nContacts(end + 1, 1) = height(subset); %#ok<AGROW>
    estimate(end + 1, 1) = stats.Estimate(index); %#ok<AGROW>
    low(end + 1, 1) = stats.Lower(index); %#ok<AGROW>
    high(end + 1, 1) = stats.Upper(index); %#ok<AGROW>
    pValue(end + 1, 1) = stats.pValue(index); %#ok<AGROW>
    sdUnits(end + 1, 1) = stats.Estimate(index) / std(subset.Coupling); %#ok<AGROW>
end

robustness = table(name, nContacts, estimate, low, high, pValue, sdUnits, ...
    'VariableNames', {'specification', 'n_contacts', 'estimate', ...
    'ci95_low', 'ci95_high', 'p_value', 'standardised'});
phg.writeTableAtomic(robustness, ...
    fullfile(cfg.tableDir, 'continuous_robustness.csv'));

fprintf(['[PHG] Continuous robustness: %s\n'], strjoin(arrayfun(@(k) ...
    sprintf('%s %.2f SD', robustness.specification(k), ...
    robustness.standardised(k)), (1:height(robustness))', ...
    'UniformOutput', false)', '; '));

outputs = struct('robustness', robustness);
end
