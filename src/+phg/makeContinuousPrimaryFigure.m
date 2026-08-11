function qc = makeContinuousPrimaryFigure(T, cfg, options)
%MAKECONTINUOUSPRIMARYFIGURE The primary result, on the scale it is stated on.
%
% The paper's primary analysis moved to a continuous signed coupling value, but
% its only figure still showed the superseded binary one: counts of
% dilation-linked contacts, per-patient dilation fractions, and sensitivity to
% a threshold the primary no longer applies. The primary result was
% unillustrated and the text pointed at a figure of a different analysis.
%
% Four panels, matching the four claims the Results now make:
%   a  the distribution of coupling by contact class, which is the primary
%   b  the region ordering, monotonic from mesial temporal to neocortex
%   c  the within-shaft spatial gradient
%   d  the frequency dependence, which is where the two structures diverge

arguments
    T table
    cfg (1,1) struct
    options.wavelet table = table()
    options.stem (1,1) string = "Fig4b_continuous_primary"
end

style = phg.figureStyle;
region = phg.cleanRegionLabels(T.NMM);
scale = median(abs(T.FitHeight), 'omitnan');
coupling = asinh(T.FitHeight ./ scale);
isHippocampal = region == "Hippocampus";
usable = isfinite(coupling);

figureHandle = figure('Units', 'inches', 'Position', [1 1 7.2 5.8], ...
    'Color', 'w');
layout = tiledlayout(figureHandle, 2, 2);
phg.setSafeLayout(layout);

%% ------------------------------------------------------------- panel a
ax = nexttile(layout);
hold(ax, 'on');
edges = linspace(-4, 4, 45);
histogram(ax, coupling(usable & ~isHippocampal), edges, ...
    'Normalization', 'probability', 'FaceColor', style.dilation, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.65);
histogram(ax, coupling(usable & isHippocampal), edges, ...
    'Normalization', 'probability', 'FaceColor', style.constriction, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.75);
yLimits = ylim(ax);
plot(ax, [0 0], yLimits, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
plot(ax, median(coupling(usable & ~isHippocampal), 'omitnan') * [1 1], ...
    yLimits, '--', 'Color', style.dilation, 'LineWidth', 1.2);
plot(ax, median(coupling(usable & isHippocampal), 'omitnan') * [1 1], ...
    yLimits, '--', 'Color', style.constriction, 'LineWidth', 1.2);
ylim(ax, yLimits);
hold(ax, 'off');
xlabel(ax, 'Coupling (asinh signed amplitude)');
ylabel(ax, 'Proportion of contacts');
title(ax, sprintf('All %d contacts, no selection', sum(usable)));
legend(ax, {sprintf('Extrahippocampal (%d)', sum(usable & ~isHippocampal)), ...
    sprintf('Hippocampal (%d)', sum(usable & isHippocampal))}, ...
    'Location', 'northwest', 'Box', 'off', 'FontSize', style.tickFontSize);
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'a');

%% ------------------------------------------------------------- panel b
ax = nexttile(layout);
regionTable = readtable(fullfile(cfg.tableDir, ...
    'continuous_region_model.csv'), 'TextType', 'string');
regionTable = regionTable(startsWith(regionTable.term, "Region_"), :);
regionTable.label = phg.shortRegionLabel( ...
    extractAfter(regionTable.term, "Region_"));
regionTable = sortrows(regionTable, 'standardised');
y = 1:height(regionTable);

% Convert each interval from model units to outcome standard deviations. The
% ratio standardised/estimate is the same divisor for both ends of the
% interval; it is computed per row rather than once because each region's
% estimate carries its own rounding, and guarded because a region whose
% estimate is essentially zero would otherwise produce an infinite bar.
usableRow = abs(regionTable.estimate) > 1e-9;
perSd = ones(height(regionTable), 1);
perSd(usableRow) = regionTable.standardised(usableRow) ./ ...
    regionTable.estimate(usableRow);
lowSd = regionTable.ci95_low .* perSd;
highSd = regionTable.ci95_high .* perSd;
% The divisor is negative for constriction-coupled regions, which swaps the
% ends of the interval; sort them back so low is left of high.
flipped = lowSd > highSd;
[lowSd(flipped), highSd(flipped)] = deal(highSd(flipped), lowSd(flipped));

hold(ax, 'on');
plot(ax, [0 0], [0.4 height(regionTable) + 0.6], '-', ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
for k = 1:height(regionTable)
    % Colour marks significance, and marks it the same way in both
    % directions. Colouring every negative estimate while requiring q < 0.05
    % of every positive one implied a constriction finding for regions whose
    % interval comfortably spans zero.
    colour = style.gray;
    if regionTable.fdr_q_value(k) < 0.05
        if regionTable.standardised(k) < 0
            colour = style.constriction;
        else
            colour = style.dilation;
        end
    end
    plot(ax, [lowSd(k) highSd(k)], [y(k) y(k)], '-', ...
        'Color', colour, 'LineWidth', 1.1);
    plot(ax, regionTable.standardised(k), y(k), 'o', 'MarkerSize', 4.5, ...
        'MarkerFaceColor', colour, 'MarkerEdgeColor', 'none');
end
hold(ax, 'off');
set(ax, 'YTick', y, 'YTickLabel', regionTable.label);
ylim(ax, [0.4, height(regionTable) + 0.6]);
% Limits follow the intervals, so that no interval is silently truncated at
% the axis edge and read as narrower than it is.
span = [min(lowSd), max(highSd)];
padding = 0.06 * diff(span);
xlim(ax, [span(1) - padding, span(2) + padding]);
xlabel(ax, 'Coupling relative to rest of implant (SD)');
title(ax, 'Region ordering on the primary scale');
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'b');

%% ------------------------------------------------------------- panel c
ax = nexttile(layout);
coordinates = T.XYZMNI;
centroid = mean(coordinates(isHippocampal & all(isfinite(coordinates), 2), :), 1);
distance = sqrt(sum((coordinates - centroid) .^ 2, 2));
lead = string(T.PtID) + "/" + phg.parseLeadLabel(T.Label);
keep = usable & all(isfinite(coordinates), 2);
[shaftId, ~] = findgroups(lead(keep));
shaftMean = splitapply(@mean, distance(keep), shaftId);
within = distance(keep) - shaftMean(shaftId);
values = coupling(keep);

binEdges = linspace(prctile(within, 2), prctile(within, 98), 9);
centres = movmean(binEdges, 2, 'Endpoints', 'discard');
binned = nan(numel(centres), 1);
spread = nan(numel(centres), 1);
for k = 1:numel(centres)
    inBin = within >= binEdges(k) & within < binEdges(k + 1);
    if sum(inBin) > 5
        binned(k) = mean(values(inBin));
        spread(k) = std(values(inBin)) / sqrt(sum(inBin));
    end
end
hold(ax, 'on');
plot(ax, [min(binEdges) max(binEdges)], [0 0], '-', ...
    'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
errorbar(ax, centres, binned, spread, 'o-', 'Color', style.constriction, ...
    'MarkerFaceColor', style.constriction, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 4.5, 'LineWidth', 1.2, 'CapSize', 3);
hold(ax, 'off');
xlabel(ax, 'Distance from hippocampus along shaft (mm, centred)');
ylabel(ax, 'Coupling');
title(ax, 'Within-shaft spatial gradient');
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'c');

%% ------------------------------------------------------------- panel d
ax = nexttile(layout);
if ~isempty(options.wavelet)
    W = options.wavelet;
    regionW = phg.cleanRegionLabels(W.NMM);
    hippW = regionW == "Hippocampus";
    scaleW = median(abs(W.FitHeight), 'omitnan');
    couplingW = asinh(W.FitHeight ./ scaleW);
    frequencies = unique(W.freq);
    hippProfile = nan(numel(frequencies), 1);
    otherProfile = nan(numel(frequencies), 1);
    for k = 1:numel(frequencies)
        atFrequency = W.freq == frequencies(k) & isfinite(couplingW);
        hippProfile(k) = median(couplingW(atFrequency & hippW), 'omitnan');
        otherProfile(k) = median(couplingW(atFrequency & ~hippW), 'omitnan');
    end
    hold(ax, 'on');
    plot(ax, frequencies, zeros(size(frequencies)), '-', ...
        'Color', [0.6 0.6 0.6], 'LineWidth', 0.6);
    plot(ax, frequencies, movmean(otherProfile, 5), '-', ...
        'Color', style.dilation, 'LineWidth', 1.4);
    plot(ax, frequencies, movmean(hippProfile, 5), '-', ...
        'Color', style.constriction, 'LineWidth', 1.4);
    hold(ax, 'off');
    set(ax, 'XScale', 'log');
    xlim(ax, [min(frequencies) max(frequencies)]);
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'Coupling');
    title(ax, 'Frequency profile by structure');
    legend(ax, {'', 'Extrahippocampal', 'Hippocampal'}, 'Location', ...
        'northwest', 'Box', 'off', 'FontSize', style.tickFontSize);
end
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'd');

qc = phg.exportPublicationFigure(figureHandle, cfg, options.stem);
close(figureHandle);
end
