function makeDerivedPublicationFigures(tables, cfg, analysis)
%MAKEDERIVEDPUBLICATIONFIGURES Rebuild manuscript figures in MATLAB only.
%
%   Figures are built around the hurdle decomposition in phg.runDerivedAnalyses.
%   "Excursion contacts" are contacts with a non-zero stored signed response,
%   i.e. contacts on which the peri-peak pupil trace left the surrogate noise
%   band at all. That is a threshold-free definition. "Legacy-selected"
%   contacts are the subset passing the historical RespSig > 0.10 contiguity
%   criterion and are only ever shown as a descriptive overlay.

arguments
    tables (1,1) struct
    cfg (1,1) struct
    analysis (1,1) struct
end

channel = tables.channel;
channel.Region = phg.cleanRegionLabels(channel.NMM);
channel.LegacySelected = channel.RespSig > cfg.legacySelectionThreshold;
channel.HasExcursion = channel.RespAreaNet ~= 0;
channel.IsDilation = channel.RespAreaNet > 0;

channelPair = tables.channelPair;
channelPair.Region = phg.cleanRegionLabels(channelPair.NMM);
channelPair.SourceRegion = phg.cleanRegionLabels(channelPair.NMM2);
channelPair.LegacySelected = channelPair.RespSig > cfg.legacySelectionThreshold;

wavelet = tables.wavelet;
wavelet.Region = phg.cleanRegionLabels(wavelet.NMM);
wavelet.LegacySelected = wavelet.RespSig > cfg.legacySelectionThreshold;

figureRegionalExtent(channel, cfg, analysis);
figurePolarity(channel, cfg, analysis);
figureSpectral(wavelet, cfg);
figureNetworkTemporal(channel, channelPair, wavelet, cfg);
figureTaskComparison(channel, cfg);
figurePeriPeak(channel, cfg);
end

function figureRegionalExtent(T, cfg, analysis)
style = phg.figureStyle;
[group, region] = findgroups(T.Region);
nChannels = splitapply(@numel, T.Chan, group);
nExcursion = splitapply(@sum, T.HasExcursion, group);
nPatients = splitapply(@(x) numel(unique(x)), string(T.PtID), group);
nExcursionPatients = splitapply(@(x,s) numel(unique(x(logical(s)))), ...
    string(T.PtID), T.HasExcursion, group);
summary = table(region, nChannels, nExcursion, nPatients, nExcursionPatients);

% Both panels must describe the same regions. Panel b can only show regions
% that entered the model as named levels; everything else was pooled into the
% "Other" reference category. Restricting panel a to that same set avoids
% showing a region in the counts panel that then has no modelled estimate.
modelledRegions = extractAfter( ...
    analysis.prevalenceTable.term(analysis.prevalenceTable.in_fdr_family), ...
    "Region_");
summary = summary(ismember(summary.region, modelledRegions), :);
summary = sortrows(summary, 'nExcursion', 'descend');

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 18.3 9.0]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
phg.setSafeLayout(layout);

axA = nexttile(layout, 1);
y = (1:height(summary))';
bars = barh(axA, y, summary.nExcursion, 0.72, 'FaceColor', 'flat', ...
    'EdgeColor', 'none');
bars.CData = phg.regionColor(summary.region);
axA.YTick = y;
axA.YTickLabel = phg.shortRegionLabel(summary.region);
axA.YDir = 'reverse';
xMaximum = max(summary.nExcursion) * 1.48;
xlim(axA, [0 xMaximum]);
ylim(axA, [0.4 height(summary)+0.6]);
for k = 1:height(summary)
    text(axA, summary.nExcursion(k) + 0.015*xMaximum, k, ...
        sprintf('%d/%d  (%d/%d pts)', summary.nExcursion(k), ...
        summary.nChannels(k), summary.nExcursionPatients(k), ...
        summary.nPatients(k)), 'FontName', style.fontName, ...
        'FontSize', 6.3, 'Color', [0.30 0.30 0.30], ...
        'VerticalAlignment', 'middle');
end
xlabel(axA, 'Contacts with a pupil excursion (n)');
title(axA, 'Regional coverage and coupling prevalence');
phg.styleAxes(axA);
phg.addPanelLabel(axA, 'a');

% Panel b: prevalence odds ratios from the hurdle part-1 model, so the reader
% sees the modelled effect rather than a raw proportion confounded by coverage.
axB = nexttile(layout, 2);
prevalence = analysis.prevalenceTable;
prevalence = prevalence(prevalence.in_fdr_family, :);
prevalence.regionName = extractAfter(prevalence.term, "Region_");
prevalence = sortrows(prevalence, 'odds_ratio', 'descend');
hold(axB, 'on');
yPosition = (1:height(prevalence))';
for k = 1:height(prevalence)
    color = phg.regionColor(prevalence.regionName(k));
    plot(axB, [prevalence.odds_ratio_ci95_low(k) ...
        prevalence.odds_ratio_ci95_high(k)], [k k], '-', ...
        'Color', color, 'LineWidth', 1.1);
    isSignificant = prevalence.fdr_q_value(k) < 0.05;
    scatter(axB, prevalence.odds_ratio(k), k, 26, color, 'o', ...
        'filled', 'MarkerEdgeColor', 'white', 'LineWidth', 0.4, ...
        'MarkerFaceAlpha', 0.35 + 0.65*isSignificant);
    if isSignificant
        text(axB, prevalence.odds_ratio_ci95_high(k)*1.18, k, ...
            sprintf('q=%.3g', prevalence.fdr_q_value(k)), ...
            'FontName', style.fontName, 'FontSize', 6.0, ...
            'Color', [0.30 0.30 0.30], 'VerticalAlignment', 'middle');
    end
end
xline(axB, 1, ':', 'Color', [0.45 0.45 0.45]);
set(axB, 'XScale', 'log');
axB.YTick = yPosition;
axB.YTickLabel = phg.shortRegionLabel(prevalence.regionName);
axB.YDir = 'reverse';
ylim(axB, [0.4 height(prevalence)+0.6]);
xlim(axB, [0.05 90]);
xlabel(axB, 'Odds of an excursion vs rest of implant');
title(axB, 'Coupling prevalence (mixed model)');
phg.styleAxes(axB);
phg.addPanelLabel(axB, 'b');

phg.exportPublicationFigure(fig, cfg, "Fig1_regional_extent_matlab");
close(fig);
end

function figurePolarity(T, cfg, analysis)
%FIGUREPOLARITY Headline figure: the hippocampal polarity reversal.
style = phg.figureStyle;
excursion = T(T.HasExcursion, :);
[group, region] = findgroups(excursion.Region);
dilation = splitapply(@(x) sum(x > 0), excursion.RespAreaNet, group);
constriction = splitapply(@(x) sum(x < 0), excursion.RespAreaNet, group);
summary = table(region, dilation, constriction);
summary.total = summary.dilation + summary.constriction;
summary = sortrows(summary, 'total', 'descend');
summary = summary(summary.total >= 5, :);
summary = summary(1:min(11,height(summary)), :);

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 18.3 7.6]);
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
phg.setSafeLayout(layout);

% ---- Panel a: composition by region
ax = nexttile(layout, 1);
y = (1:height(summary))';
bars = barh(ax, y, [summary.dilation summary.constriction], ...
    0.74, 'stacked', 'EdgeColor', 'white', 'LineWidth', 0.25);
bars(1).FaceColor = style.dilation;
bars(2).FaceColor = style.constriction;
ax.YTick = y;
ax.YTickLabel = phg.shortRegionLabel(summary.region);
ax.YDir = 'reverse';
ylim(ax, [0.4 height(summary)+0.6]);
xlim(ax, [0 max(summary.total)*1.06]);
for k = 1:height(summary)
    if summary.dilation(k) >= 3
        text(ax, summary.dilation(k)/2, k, string(summary.dilation(k)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontName', style.fontName, 'FontSize', 6.2, ...
            'FontWeight', 'bold', 'Color', 'white');
    end
    if summary.constriction(k) >= 3
        text(ax, summary.dilation(k)+summary.constriction(k)/2, k, ...
            string(summary.constriction(k)), 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontName', style.fontName, ...
            'FontSize', 6.2, 'FontWeight', 'bold', 'Color', 'white');
    end
end
xlabel(ax, 'Excursion contacts (n)');
title(ax, 'Coupling direction by region');
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'a');
legend(ax, bars, {'Dilation', 'Constriction'}, 'Location', 'southeast', ...
    'FontSize', 6.5, 'Box', 'off');

% ---- Panel b: within-patient paired rates
axB = nexttile(layout, 2);
paired = analysis.pairedTable;
paired = paired(~isnan(paired.hippocampal_dilation_rate) & ...
    ~isnan(paired.extrahippocampal_dilation_rate), :);
% Panel a's legend defines coral as dilation-linked and teal as
% constriction-linked. Re-using those two colours here for contact class and
% for the direction of a patient's change gave one colour three meanings in a
% single figure, and put a coral -- "dilation" -- marker at a dilation
% fraction of zero. Both axes of this panel already encode contact class by
% position, so the markers carry no colour information and the only colour
% left is the one distinction a reader needs: which patients run against the
% predicted direction.
hold(axB, 'on');
for k = 1:height(paired)
    difference = paired.hippocampal_dilation_rate(k) - ...
        paired.extrahippocampal_dilation_rate(k);
    if difference < 0
        lineColor = style.gray;
        lineWidth = 0.9;
    else
        lineColor = style.sensorimotor;
        lineWidth = 1.5;
    end
    plot(axB, [1 2], [paired.extrahippocampal_dilation_rate(k) ...
        paired.hippocampal_dilation_rate(k)], '-', 'Color', ...
        [lineColor 0.75], 'LineWidth', lineWidth);
end
scatter(axB, ones(height(paired),1), ...
    paired.extrahippocampal_dilation_rate, 24, style.gray, 'filled', ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.4, 'MarkerFaceAlpha', 0.9);
scatter(axB, 2*ones(height(paired),1), ...
    paired.hippocampal_dilation_rate, 24, style.gray, 'filled', ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.4, 'MarkerFaceAlpha', 0.9);
axB.XTick = [1 2];
axB.XTickLabel = {'Extra-hipp.', 'Hippocampus'};
xlim(axB, [0.65 2.35]);
ylim(axB, [-0.06 1.18]);
ylabel(axB, 'Dilation-linked fraction');
title(axB, 'Within-patient comparison');
lowerCount = sum(paired.dilation_rate_difference < 0);
informative = sum(paired.dilation_rate_difference ~= 0);
signTestP = analysis.designSummary.value( ...
    analysis.designSummary.parameter == "paired_sign_test_p_value");
% Wrapped onto two lines so the annotation stays inside its own tile instead of
% spilling into the neighbouring panel.
text(axB, 1.5, 1.14, {sprintf('%d/%d patients lower', lowerCount, informative), ...
    sprintf('sign test p = %.3f', signTestP)}, 'Units', 'data', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontName', style.fontName, 'FontSize', 6.0, 'Color', [0.30 0.30 0.30]);
phg.styleAxes(axB);
phg.addPanelLabel(axB, 'b');

% ---- Panel c: threshold sensitivity
axC = nexttile(layout, 3);
sensitivity = analysis.polaritySensitivity;
valid = isfinite(sensitivity.odds_ratio_hippocampal);
sensitivity = sensitivity(valid, :);
hold(axC, 'on');
for k = 1:height(sensitivity)
    plot(axC, [sensitivity.odds_ratio_ci95_low(k) ...
        sensitivity.odds_ratio_ci95_high(k)], [k k], '-', ...
        'Color', style.constriction, 'LineWidth', 1.1);
end
scatter(axC, sensitivity.odds_ratio_hippocampal, ...
    (1:height(sensitivity))', 26, style.constriction, 'filled', ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.4);
xline(axC, 1, ':', 'Color', [0.45 0.45 0.45]);
set(axC, 'XScale', 'log');
axC.YTick = 1:height(sensitivity);
axC.YTickLabel = compose('%.2f', sensitivity.selection_threshold);
axC.YDir = 'reverse';
ylim(axC, [0.4 height(sensitivity)+0.6]);
ylabel(axC, 'Contiguity selection threshold');
xlabel(axC, 'Odds of dilation, hippocampal vs not');
title(axC, 'Threshold sensitivity');
phg.styleAxes(axC);
phg.addPanelLabel(axC, 'c');

phg.exportPublicationFigure(fig, cfg, "Fig2_polarity_matlab");
close(fig);
end

function figureSpectral(W, cfg)
selected = W(W.LegacySelected & isfinite(W.RespAreaAbs) & ...
    W.RespAreaAbs > 0 & isfinite(W.freq), :);
[group, region] = findgroups(selected.Region);
nChannels = splitapply(@(p,c) numel(unique(p+"_"+string(c))), ...
    string(selected.PtID), selected.Chan, group);
nPatients = splitapply(@(x) numel(unique(x)), string(selected.PtID), group);
coverage = table(region, nChannels, nPatients);
coverage = coverage(coverage.nPatients >= 2, :);
coverage = sortrows(coverage, 'nChannels', 'descend');
keep = coverage.region(1:min(8,height(coverage)));
profileRegions = keep(1:min(5,numel(keep)));

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 18.3 13.0]);
layout = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
phg.setSafeLayout(layout);
axA = nexttile(layout, 1);
hold(axA, 'on');
lineHandles = gobjects(numel(profileRegions),1);
for r = 1:numel(profileRegions)
    sub = selected(selected.Region == profileRegions(r), :);
    [frequencyGroup, frequency] = findgroups(sub.freq);
    medianAmplitude = splitapply(@median, sub.RespAreaAbs, frequencyGroup);
    lower = splitapply(@(x) prctile(x,25), sub.RespAreaAbs, frequencyGroup);
    upper = splitapply(@(x) prctile(x,75), sub.RespAreaAbs, frequencyGroup);
    [frequency, order] = sort(frequency);
    medianAmplitude = medianAmplitude(order);
    lower = lower(order);
    upper = upper(order);
    color = phg.regionColor(profileRegions(r));
    fill(axA, [frequency; flipud(frequency)], ...
        [lower; flipud(upper)], color, 'FaceAlpha', 0.16, ...
        'EdgeColor', 'none');
    lineHandles(r) = loglog(axA, frequency, medianAmplitude, ...
        'Color', color, 'LineWidth', 1.35);
end
set(axA, 'XScale', 'log', 'YScale', 'log');
xlim(axA, [2 256]);
xlabel(axA, 'Frequency (Hz)');
ylabel(axA, 'Pupil-coupling amplitude (a.u.)');
title(axA, 'Spectral profiles (median and interquartile range)');
phg.styleAxes(axA);
phg.addPanelLabel(axA, 'a');
legendHandle = legend(axA, lineHandles, ...
    cellstr(phg.shortRegionLabel(profileRegions)), ...
    'Orientation', 'horizontal', 'NumColumns', numel(profileRegions), ...
    'FontSize', 7, 'Box', 'off');
legendHandle.Layout.Tile = 'south';

axB = nexttile(layout, 2);
slope = nan(numel(keep),1);
fitR2 = nan(numel(keep),1);
for r = 1:numel(keep)
    sub = selected(selected.Region == keep(r), :);
    [frequencyGroup, frequency] = findgroups(sub.freq);
    amplitude = splitapply(@median, sub.RespAreaAbs, frequencyGroup);
    valid = frequency >= 2 & frequency <= 256 & amplitude > 0 & ...
        isfinite(frequency) & isfinite(amplitude);
    x = log10(frequency(valid));
    y = log10(amplitude(valid));
    coefficients = robustfit(x, y);
    estimate = coefficients(1) + coefficients(2).*x;
    slope(r) = coefficients(2);
    fitR2(r) = corr(y, estimate, 'Rows', 'complete').^2;
end
slopeTable = table(keep, slope, fitR2, ...
    'VariableNames', {'region', 'robust_loglog_slope', 'fit_r_squared'});
slopeTable = sortrows(slopeTable, 'robust_loglog_slope');
yPosition = (1:height(slopeTable))';
bars = barh(axB, yPosition, slopeTable.robust_loglog_slope, ...
    0.72, 'FaceColor', 'flat', 'EdgeColor', 'none');
bars.CData = phg.regionColor(slopeTable.region);
axB.YTick = yPosition;
axB.YTickLabel = phg.shortRegionLabel(slopeTable.region);
xline(axB, 0, ':', 'Color', [0.45 0.45 0.45]);
xlabel(axB, 'Robust log-log spectral slope');
title(axB, 'Robust log-log trend summary');
phg.styleAxes(axB);
phg.addPanelLabel(axB, 'b');
phg.writeTableAtomic(slopeTable, fullfile(cfg.tableDir, ...
    'spectral_robust_slopes.csv'));

phg.exportPublicationFigure(fig, cfg, "Fig3_spectral_profiles_matlab");
close(fig);
end

function figureNetworkTemporal(T, P, W, cfg)
style = phg.figureStyle;
selectedPairs = P(P.LegacySelected, :);
[group, source, target] = findgroups(selectedPairs.SourceRegion, ...
    selectedPairs.Region);
weight = splitapply(@numel, selectedPairs.Chan, group);
nPatients = splitapply(@(x) numel(unique(x)), ...
    string(selectedPairs.PtID), group);
edgeTable = table(source, target, weight, nPatients);
edgeTable = edgeTable(edgeTable.weight >= 3 & edgeTable.nPatients >= 2 & ...
    edgeTable.source ~= edgeTable.target, :);
G = digraph(cellstr(edgeTable.source), cellstr(edgeTable.target), ...
    edgeTable.weight);
nodeNames = string(G.Nodes.Name);
pageRank = centrality(G, 'pagerank', 'Importance', G.Edges.Weight);
centralityTable = table(nodeNames, pageRank, ...
    'VariableNames', {'region', 'pagerank'});
centralityTable = sortrows(centralityTable, 'pagerank', 'descend');
phg.writeTableAtomic(edgeTable, fullfile(cfg.tableDir, ...
    'high_gamma_network_edges.csv'));
phg.writeTableAtomic(centralityTable, fullfile(cfg.tableDir, ...
    'high_gamma_network_centrality.csv'));

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 18.3 11.5]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
phg.setSafeLayout(layout);
axA = nexttile(layout, 1, [2 1]);
plotNetworkPanel(axA, G, centralityTable);
title(axA, 'Descriptive high-gamma pair network');
phg.addPanelLabel(axA, 'a');

axB = nexttile(layout, 2);
lagGroups = {
    "MTG middle temporal gyrus", "MTG", phg.regionColor("MTG middle temporal gyrus")
    "Hippocampus", "Hipp", phg.regionColor("Hippocampus")
    ["PoG postcentral gyrus","PrG precentral gyrus"], "S/M", style.sensorimotor
    };
hold(axB, 'on');
allLag = zeros(0,1);
allGroup = strings(0,1);
for r = 1:size(lagGroups,1)
    mask = T.LegacySelected & ismember(T.Region, lagGroups{r,1}) & ...
        T.FitR2 > 0.5 & abs(T.FitShift) < 5000;
    lag = T.FitShift(mask) ./ 1000;
    allLag = [allLag; lag]; %#ok<AGROW>
    allGroup = [allGroup; repmat(lagGroups{r,2}, numel(lag), 1)]; %#ok<AGROW>
    rng(300+r, 'twister');
    jitter = (rand(size(lag))-0.5).*0.24;
    scatter(axB, r+jitter, lag, 15, lagGroups{r,3}, 'filled', ...
        'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'white', ...
        'LineWidth', 0.25);
    q = prctile(lag, [25 50 75]);
    plot(axB, [r-0.20 r+0.20], [q(2) q(2)], 'Color', lagGroups{r,3}, ...
        'LineWidth', 1.7);
    plot(axB, [r r], [q(1) q(3)], 'Color', lagGroups{r,3}, ...
        'LineWidth', 1.0);
end
yline(axB, 0, ':', 'Color', [0.50 0.50 0.50]);
axB.XTick = 1:3;
axB.XTickLabel = lagGroups(:,2);
axB.YLim = [-3 4];
ylabel(axB, 'Pupil lag from HG peak (s)');
title(axB, 'Stored fit-lag distribution');
phg.styleAxes(axB);
phg.addPanelLabel(axB, 'b');
try
    [pValue, analysisTable] = kruskalwallis(allLag, cellstr(allGroup), 'off');
    hStatistic = analysisTable{2,5};
    text(axB, 0.98, 0.04, sprintf('Kruskal H=%.1f, p=%.2g', ...
        hStatistic, pValue), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'FontSize', 6.3, ...
        'Color', [0.30 0.30 0.30]);
catch
end

axC = nexttile(layout, 4);
selectedWavelet = W(W.LegacySelected & W.FitR2 > 0.5 & ...
    abs(W.FitShift) < 5000, :);
focus = ["MTG middle temporal gyrus", "Hippocampus", ...
    "PrG precentral gyrus", "PoG postcentral gyrus"];
hold(axC, 'on');
lineHandles = gobjects(numel(focus),1);
for r = 1:numel(focus)
    sub = selectedWavelet(selectedWavelet.Region == focus(r), :);
    [frequencyGroup, frequency] = findgroups(sub.freq);
    lag = splitapply(@median, sub.FitShift./1000, frequencyGroup);
    [frequency, order] = sort(frequency);
    lag = lag(order);
    lineHandles(r) = semilogx(axC, frequency, lag, ...
        'Color', phg.regionColor(focus(r)), 'LineWidth', 1.2);
end
yline(axC, 0, ':', 'Color', [0.50 0.50 0.50]);
xlim(axC, [2 256]);
ylim(axC, [-5 5]);
xlabel(axC, 'Frequency (Hz)');
ylabel(axC, 'Median lag (s)');
title(axC, 'Frequency-resolved stored fit lag');
phg.styleAxes(axC);
phg.addPanelLabel(axC, 'c');
legend(axC, lineHandles, cellstr(phg.shortRegionLabel(focus)), ...
    'Location', 'northeast', 'Orientation', 'horizontal', ...
    'NumColumns', 2, 'FontSize', 6.5, 'Box', 'off');

phg.exportPublicationFigure(fig, cfg, "Fig4_network_temporal_matlab");
close(fig);
end

function plotNetworkPanel(ax, G, centralityTable)
nodeNames = string(G.Nodes.Name);
hubMTG = "MTG middle temporal gyrus";
hubHipp = "Hippocampus";
ends = string(G.Edges.EndNodes);
weights = G.Edges.Weight;
candidate = strings(0,1);
cluster = strings(0,1);
strength = zeros(0,1);
for k = 1:numel(nodeNames)
    if ismember(nodeNames(k), [hubMTG hubHipp])
        continue
    end
    mtgWeight = sum(weights((ends(:,1)==nodeNames(k) & ends(:,2)==hubMTG) | ...
        (ends(:,2)==nodeNames(k) & ends(:,1)==hubMTG)));
    hippWeight = sum(weights((ends(:,1)==nodeNames(k) & ends(:,2)==hubHipp) | ...
        (ends(:,2)==nodeNames(k) & ends(:,1)==hubHipp)));
    if max(mtgWeight, hippWeight) > 0
        candidate(end+1,1) = nodeNames(k); %#ok<AGROW>
        strength(end+1,1) = max(mtgWeight, hippWeight); %#ok<AGROW>
        if mtgWeight >= hippWeight
            cluster(end+1,1) = "MTG"; %#ok<AGROW>
        else
            cluster(end+1,1) = "Hipp"; %#ok<AGROW>
        end
    end
end
keep = [hubMTG; hubHipp];
for name = ["MTG" "Hipp"]
    index = find(cluster == name);
    [~, order] = sort(strength(index), 'descend');
    index = index(order(1:min(4,numel(order))));
    keep = [keep; candidate(index)]; %#ok<AGROW>
end
keep = unique(keep, 'stable');
keep = keep(ismember(keep, nodeNames));
Gp = subgraph(G, cellstr(keep));
names = string(Gp.Nodes.Name);
x = zeros(numel(names),1);
y = zeros(numel(names),1);
mtgSatellite = names(~ismember(names,[hubMTG hubHipp]) & ...
    ismember(names, candidate(cluster=="MTG")));
hippSatellite = names(~ismember(names,[hubMTG hubHipp]) & ...
    ismember(names, candidate(cluster=="Hipp")));
x(names==hubMTG) = -0.72;
x(names==hubHipp) = 0.72;
theta = linspace(2.2, 4.1, max(1,numel(mtgSatellite)));
for k = 1:numel(mtgSatellite)
    index = names == mtgSatellite(k);
    x(index) = -0.72 + 0.78*cos(theta(k));
    y(index) = 0.95*sin(theta(k));
end
theta = linspace(0.95, -0.95, max(1,numel(hippSatellite)));
for k = 1:numel(hippSatellite)
    index = names == hippSatellite(k);
    x(index) = 0.72 + 0.78*cos(theta(k));
    y(index) = 0.95*sin(theta(k));
end
nodeCentrality = zeros(numel(names),1);
for k = 1:numel(names)
    nodeCentrality(k) = centralityTable.pagerank( ...
        centralityTable.region == names(k));
end
edgeWidth = 0.5 + 2.5 .* Gp.Edges.Weight ./ max(Gp.Edges.Weight);
graphHandle = plot(ax, Gp, 'XData', x, 'YData', y, ...
    'NodeLabel', repmat({''}, numel(names), 1), ...
    'NodeColor', phg.regionColor(names), ...
    'MarkerSize', 5 + 26.*nodeCentrality./max(nodeCentrality), ...
    'EdgeColor', [0.55 0.55 0.55], 'LineWidth', edgeWidth, ...
    'ArrowSize', 8);
graphHandle.EdgeAlpha = 0.46;
hold(ax, 'on');
for k = 1:numel(names)
    if ismember(names(k), [hubMTG hubHipp])
        text(ax, x(k), y(k), phg.shortRegionLabel(names(k)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'Color', 'white', 'FontWeight', 'bold', 'FontSize', 7);
    else
        direction = sign(x(k));
        if direction == 0, direction = 1; end
        text(ax, x(k)+0.08*direction, y(k), phg.shortRegionLabel(names(k)), ...
            'HorizontalAlignment', ternary(direction<0,'right','left'), ...
            'VerticalAlignment', 'middle', 'FontSize', 6.3);
    end
end
axis(ax, 'normal', 'off');
% The node coordinates span roughly [-1.9, 1.9] once satellite labels are
% placed, so the limits are set from the drawn extent rather than padded to a
% fixed box that would leave a blank column on each side.
xlim(ax, [min(x)-0.42 max(x)+0.42]);
ylim(ax, [min(y)-0.12 max(y)+0.12]);
phg.styleAxes(ax);
end

function value = ternary(condition, trueValue, falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end

function figureTaskComparison(T, cfg)
style = phg.figureStyle;
task = taskForPatient(string(T.PtID), cfg);
[group, patient, taskGroup] = findgroups(string(T.PtID), task);
nContacts = splitapply(@numel, T.Chan, group);
nExcursion = splitapply(@sum, T.HasExcursion, group);
patientTable = table(patient, taskGroup, nContacts, nExcursion, ...
    nExcursion./nContacts, 'VariableNames', {'patient', 'task', ...
    'n_contacts', 'n_excursion_contacts', 'proportion'});
order = ["Working memory", "Thermal perception", "Rest"];

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 8.9 8.4]);
layout = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
phg.setSafeLayout(layout);
ax = nexttile(layout);
hold(ax, 'on');
allValue = zeros(0,1);
allGroup = strings(0,1);
tickLabel = strings(3,1);
for k = 1:numel(order)
    values = patientTable.proportion(patientTable.task == order(k));
    allValue = [allValue; values]; %#ok<AGROW>
    allGroup = [allGroup; repmat(order(k),numel(values),1)]; %#ok<AGROW>
    % The interquartile box is drawn first. Drawn afterwards its opaque face
    % hides every patient whose proportion falls inside the box, which on this
    % panel was two thirds of the working-memory group.
    q = prctile(values, [25 50 75]);
    lightColor = 0.84.*[1 1 1] + 0.16.*style.task(k,:);
    rectangle(ax, 'Position', [k-0.25 q(1) 0.50 q(3)-q(1)], ...
        'FaceColor', lightColor, ...
        'EdgeColor', style.task(k,:), 'LineWidth', 0.8);
    plot(ax, [k-0.25 k+0.25], [q(2) q(2)], ...
        'Color', style.task(k,:), 'LineWidth', 1.5);
    rng(500+k, 'twister');
    jitter = (rand(size(values))-0.5).*0.22;
    scatter(ax, k+jitter, values, 28, style.task(k,:), 'filled', ...
        'MarkerFaceAlpha', 0.82, 'MarkerEdgeColor', 'white', ...
        'LineWidth', 0.35);
    tickLabel(k) = order(k);
    text(ax, k, max(values)+0.012, sprintf('n=%d', numel(values)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontName', style.fontName, 'FontSize', 6.3, ...
        'Color', [0.35 0.35 0.35]);
end
ax.XTick = 1:3;
ax.XTickLabel = tickLabel;
ax.XLim = [0.5 3.5];
ax.YLim = [0 max(0.45, max(allValue)*1.18)];
ylabel(ax, 'Per-patient excursion proportion');
ax.XTickLabelRotation = 15;
title(ax, 'Coupling prevalence by recording condition');
phg.styleAxes(ax);
% The test annotated on the panel is written out as well, so that the
% manuscript can quote it without any number being retyped by hand.
[pValue, analysisTable] = kruskalwallis(allValue, cellstr(allGroup), 'off');
hStatistic = analysisTable{2,5};
text(ax, 0.98, 0.96, sprintf('Kruskal H=%.2f, p=%.2g', ...
    hStatistic, pValue), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 6.5, 'Color', [0.30 0.30 0.30]);
phg.writeTableAtomic(table("Kruskal-Wallis", hStatistic, ...
    numel(unique(allGroup))-1, numel(allValue), pValue, ...
    'VariableNames', {'test','statistic','df','n_patients','p_value'}), ...
    fullfile(cfg.tableDir, 'task_condition_test.csv'));
phg.writeTableAtomic(patientTable, fullfile(cfg.tableDir, ...
    'task_patient_selection_summary.csv'));
phg.exportPublicationFigure(fig, cfg, "Fig5_task_comparison_matlab");
close(fig);
end

function task = taskForPatient(patient, cfg)
%TASKFORPATIENT Look up each participant's recording condition.
%
%   The mapping lives in config/patient_task_map.tsv rather than in this file.
%   Which participant performed which task is participant-level information,
%   and the public mirror of this repository must carry none; keeping it in a
%   data file means the mirror can simply omit that file.

mapFile = fullfile(cfg.repoRoot, 'config', 'patient_task_map.tsv');
if ~isfile(mapFile)
    error('phg:missingTaskMap', ...
        ['%s is missing. It maps each participant to a recording ' ...
         'condition and is not distributed with the public code mirror.'], ...
        mapFile);
end
map = readtable(mapFile, 'FileType', 'text', 'Delimiter', '\t', ...
    'TextType', 'string');
task = repmat("Unknown", size(patient));
[found, where] = ismember(patient, map.patient_id);
task(found) = map.task(where(found));
end

function figurePeriPeak(T, cfg)
style = phg.figureStyle;
time = linspace(-20,20,2001);
timeMs = time .* 1000;
halfWindowMs = 20000;

% A Gaussian whose fitted width exceeds the analysis window is not resolvable
% from a constant within that window: it renders as a flat line at its
% amplitude and contributes a pure offset to the group mean. One MTG contact
% has a fitted width of 192 s with an amplitude of 79, which alone shifted the
% 48-contact MTG mean by roughly 1.6 units and compressed every other panel.
% Such fits are excluded (4 of 160 selected contacts) rather than displayed.
resolvableFit = T.LegacySelected & T.FitR2 > 0.5 & ...
    isfinite(T.FitWidth) & abs(T.FitWidth) > 0 & ...
    abs(T.FitWidth) <= halfWindowMs & abs(T.FitShift) <= halfWindowMs;
selected = T(resolvableFit, :);
groups = {
    "MTG middle temporal gyrus", "MTG", phg.regionColor("MTG middle temporal gyrus")
    "Hippocampus", "Hipp", phg.regionColor("Hippocampus")
    ["PoG postcentral gyrus","PrG precentral gyrus"], "S/M", style.sensorimotor
    };
meanTrace = cell(3,1);
semTrace = cell(3,1);
count = zeros(3,1);
medianLag = zeros(3,1);
for g = 1:3
    sub = selected(ismember(selected.Region,groups{g,1}), :);
    traces = nan(height(sub),numel(timeMs));
    for k = 1:height(sub)
        if isfinite(sub.FitWidth(k)) && sub.FitWidth(k) ~= 0
            traces(k,:) = sub.FitHeight(k) .* exp(-((timeMs - ...
                sub.FitShift(k))./sub.FitWidth(k)).^2);
        end
    end
    meanTrace{g} = mean(traces,1,'omitnan');
    semTrace{g} = std(traces,0,1,'omitnan') ./ sqrt(sum(any(isfinite(traces),2)));
    count(g) = height(sub);
    medianLag(g) = median(sub.FitShift,'omitnan')./1000;
end
allLower = cellfun(@(m,s) min(m-s,[],'omitnan'), meanTrace, semTrace);
allUpper = cellfun(@(m,s) max(m+s,[],'omitnan'), meanTrace, semTrace);
yLimits = [min(allLower) max(allUpper)];
padding = max(diff(yLimits)*0.10, 0.05);
yLimits = yLimits + [-padding padding];

fig = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 18.3 6.8]);
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
phg.setSafeLayout(layout);
for g = 1:3
    ax = nexttile(layout, g);
    hold(ax, 'on');
    color = groups{g,3};
    fill(ax, [time fliplr(time)], ...
        [meanTrace{g}-semTrace{g} fliplr(meanTrace{g}+semTrace{g})], ...
        color, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
    plot(ax, time, meanTrace{g}, 'Color', color, 'LineWidth', 1.45);
    yline(ax, 0, ':', 'Color', [0.55 0.55 0.55]);
    xline(ax, 0, ':', 'Color', [0.55 0.55 0.55]);
    xlim(ax, [-20 20]);
    ylim(ax, yLimits);
    xlabel(ax, 'Time from HG peak (s)');
    if g == 1
        ylabel(ax, 'Reconstructed pupil change (a.u.)');
    else
        ax.YTickLabel = {};
    end
    title(ax, sprintf('%s (n=%d)', groups{g,2}, count(g)), ...
        'Color', color, 'FontWeight', 'bold');
    text(ax, 0.04, 0.94, sprintf('Median lag %+.2f s', medianLag(g)), ...
        'Units', 'normalized', 'FontSize', 6.5, ...
        'Color', [0.30 0.30 0.30], 'VerticalAlignment', 'top');
    if g == 3
        text(ax, 0.98, 0.04, 'Stored Gaussian-fit parameters', ...
            'Units', 'normalized', 'HorizontalAlignment', 'right', ...
            'FontSize', 6.0, 'Color', [0.40 0.40 0.40]);
    end
    phg.styleAxes(ax);
    phg.addPanelLabel(ax, char('a'+g-1));
end
phg.exportPublicationFigure(fig, cfg, "Fig6_peri_peak_fits_matlab");
close(fig);
end
