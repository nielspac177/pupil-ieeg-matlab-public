function outputs = runReferenceClassChecks(channel, cfg)
%RUNREFERENCECLASSCHECKS Is the contrast hippocampal, or something coarser?
%
%   The primary comparison sets the hippocampus against everything else, and
%   "everything else" is a heterogeneous class. Two objections follow, and both
%   are answerable with the data in hand.
%
%   The first is that the contrast is really driven by middle temporal gyrus,
%   which contributes the most excursion contacts and is 91% dilation-linked.
%   The reference class is therefore rebuilt with MTG removed, with amygdala
%   removed, with both removed, and with precentral gyrus removed.
%
%   The second is that the effect may be mesial temporal rather than
%   hippocampal — the amygdala sits millimetres away, is implanted for the same
%   clinical reason, and is pooled into the reference class by the primary
%   model. A three-level model puts it in its own category so its position can
%   be read rather than assumed.
%
%   Both are reported as robustness, not as new primary results.

arguments
    channel table
    cfg (1,1) struct
end

region = phg.cleanRegionLabels(channel.NMM);
excursion = channel.RespAreaNet ~= 0;

frame = table;
frame.PtID = removecats(categorical(string(channel.PtID(excursion))));
frame.Lead = categorical(phg.parseLeadLabel(string(channel.Label(excursion))));
frame.IsDilation = double(channel.RespAreaNet(excursion) > 0);
frame.Region = region(excursion);

formula = 'IsDilation ~ ContactClass + (1|PtID) + (1|PtID:Lead)';

%% Reference-class robustness
names = ["all extrahippocampal (primary)", "excluding MTG", ...
    "excluding amygdala", "excluding MTG and amygdala", "excluding precentral"];
dropped = {string.empty, "MTG middle temporal gyrus", "Amygdala", ...
    ["MTG middle temporal gyrus", "Amygdala"], "PrG precentral gyrus"};

rows = cell(0, 1);
for v = 1:numel(names)
    subset = frame(frame.Region == "Hippocampus" | ...
        ~ismember(frame.Region, ["Hippocampus", dropped{v}]), :);
    subset.PtID = removecats(subset.PtID);
    class = repmat("Extrahippocampal", height(subset), 1);
    class(subset.Region == "Hippocampus") = "Hippocampal";
    subset.ContactClass = reordercats(categorical(class), ...
        {'Extrahippocampal', 'Hippocampal'});
    rows{end+1, 1} = localFit(subset, formula, names(v)); %#ok<AGROW>
end
referenceTable = vertcat(rows{:});
phg.writeTableAtomic(referenceTable, ...
    fullfile(cfg.tableDir, 'polarity_reference_class.csv'));

%% Three-level model: hippocampus, amygdala, neocortex
class = repmat("Neocortical", height(frame), 1);
class(frame.Region == "Hippocampus") = "Hippocampal";
class(frame.Region == "Amygdala") = "Amygdalar";
threeLevel = frame;
threeLevel.ContactClass = reordercats(categorical(class), ...
    {'Neocortical', 'Hippocampal', 'Amygdalar'});
model = fitglme(threeLevel, formula, 'Distribution', 'Binomial', ...
    'Link', 'logit', 'FitMethod', 'Laplace');
coefficients = model.Coefficients;

rows = cell(0, 1);
for k = 2:height(coefficients)
    rows{end+1, 1} = table(string(coefficients.Name(k)), ...
        exp(coefficients.Estimate(k)), exp(coefficients.Lower(k)), ...
        exp(coefficients.Upper(k)), coefficients.pValue(k), ...
        'VariableNames', {'term', 'odds_ratio', 'odds_ratio_ci95_low', ...
        'odds_ratio_ci95_high', 'p_value'}); %#ok<AGROW>
end
threeLevelTable = vertcat(rows{:});
phg.writeTableAtomic(threeLevelTable, ...
    fullfile(cfg.tableDir, 'polarity_three_level.csv'));

outputs = struct('referenceTable', referenceTable, ...
    'threeLevelTable', threeLevelTable);

for k = 1:height(referenceTable)
    fprintf('[PHG] Reference class %-30s OR %.3f [%.3f, %.3f], P = %.3g\n', ...
        referenceTable.reference_class(k), referenceTable.odds_ratio(k), ...
        referenceTable.odds_ratio_ci95_low(k), ...
        referenceTable.odds_ratio_ci95_high(k), referenceTable.p_value(k));
end
for k = 1:height(threeLevelTable)
    fprintf('[PHG] Three-level %-28s OR %.3f [%.3f, %.3f], P = %.3g\n', ...
        threeLevelTable.term(k), threeLevelTable.odds_ratio(k), ...
        threeLevelTable.odds_ratio_ci95_low(k), ...
        threeLevelTable.odds_ratio_ci95_high(k), threeLevelTable.p_value(k));
end
end

% -------------------------------------------------------------------------
function row = localFit(subset, formula, label)
model = fitglme(subset, formula, 'Distribution', 'Binomial', 'Link', 'logit', ...
    'FitMethod', 'Laplace');
coefficients = model.Coefficients;
k = find(startsWith(string(coefficients.Name), "ContactClass"), 1);
row = table(label, height(subset), exp(coefficients.Estimate(k)), ...
    exp(coefficients.Lower(k)), exp(coefficients.Upper(k)), ...
    coefficients.pValue(k), 'VariableNames', {'reference_class', ...
    'n_contacts', 'odds_ratio', 'odds_ratio_ci95_low', ...
    'odds_ratio_ci95_high', 'p_value'});
end
