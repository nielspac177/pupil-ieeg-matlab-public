function qc = makeMethodsFigure(cfg)
%MAKEMETHODSFIGURE Schematic of the peri-peak measurement and the hurdle model.
%
%   This figure is drawn entirely from synthetic illustrative traces. It uses
%   no participant data and depends on no result table, so it runs in a clone
%   that has neither, and it is the one figure that can be released publicly.
%   Its purpose is to make two things visible that prose states poorly:
%
%     - what the stored `RespSig` statistic actually measures, namely the
%       fraction of the peri-peak window spanned by the LONGEST CONTIGUOUS run
%       above the surrogate band, which is a contiguity criterion and not a
%       p-value; and
%     - why the resulting outcome cannot be modelled with a single regression,
%       because a contact with no suprathreshold run contributes an exact zero
%       rather than a small number.
%
%   Panel (c) is a schematic distribution, not the observed one; the observed
%   proportions are reported in the Results.

arguments
    cfg (1,1) struct
end

style = phg.figureStyle;
rng(20260806, 'twister');

fig = figure('Color', 'w', 'Units', 'centimeters', 'Visible', 'off', ...
    'Position', [2 2 18.0 11.6]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
phg.setSafeLayout(layout);

%% (a) High-gamma peak detection on a synthetic trace
axA = nexttile(layout);
hold(axA, 'on');
t = linspace(0, 30, 3000);
background = filterSmooth(randn(size(t)), 45);
background = 0.55 .* background ./ std(background);
highGamma = background + 2.55.*exp(-((t-7.4).^2)/0.35) ...
    + 2.95.*exp(-((t-15.9).^2)/0.30) + 2.25.*exp(-((t-24.2).^2)/0.40);
threshold = 1.8;
plot(axA, t, highGamma, 'Color', [0.35 0.38 0.42], 'LineWidth', 0.7);
yline(axA, threshold, '--', 'Color', style.gray, 'LineWidth', 0.8);
[peakValue, peakIdx] = findpeaks(highGamma, 'MinPeakHeight', threshold, ...
    'MinPeakDistance', 300);
plot(axA, t(peakIdx), peakValue, 'v', 'MarkerSize', 4.5, ...
    'MarkerFaceColor', style.dilation, 'MarkerEdgeColor', 'none');
text(axA, 0.35, threshold+0.12, 'detection threshold', ...
    'FontName', style.fontName, 'FontSize', 6.2, 'Color', style.gray, ...
    'VerticalAlignment', 'bottom');
xlabel(axA, 'Time (s)');
ylabel(axA, 'High-gamma envelope (z)');
title(axA, 'Detect high-gamma peaks in each contact');
axA.XLim = [0 30];
axA.YLim = [min(highGamma)-0.4, max(highGamma)+0.9];
phg.styleAxes(axA);
phg.addPanelLabel(axA, 'a');

%% (b) Peri-peak pupil response and the contiguity criterion
axB = nexttile(layout);
hold(axB, 'on');
lag = linspace(-5, 5, 1001);
response = 0.78.*exp(-((lag-0.9).^2)/1.9) - 0.10;
band = 0.30 + 0.05.*cos(lag./2);
fill(axB, [lag fliplr(lag)], [band fliplr(-band)], style.lightGray, ...
    'FaceAlpha', 0.55, 'EdgeColor', 'none');
plot(axB, lag, response, 'Color', style.dilation, 'LineWidth', 1.4);
yline(axB, 0, '-', 'Color', [0.75 0.76 0.78], 'LineWidth', 0.6);

above = response > band;
[runStart, runStop] = longestRun(above);
runLag = lag(runStart:runStop);
fill(axB, [runLag fliplr(runLag)], ...
    [response(runStart:runStop) fliplr(band(runStart:runStop))], ...
    style.dilation, 'FaceAlpha', 0.30, 'EdgeColor', 'none');
plot(axB, [lag(runStart) lag(runStop)], [-0.42 -0.42], '-', ...
    'Color', style.dilation, 'LineWidth', 2.0);
text(axB, mean(runLag), -0.47, ...
    sprintf('longest contiguous run = %.0f%% of window', ...
    100*(runStop-runStart+1)/numel(lag)), ...
    'FontName', style.fontName, 'FontSize', 6.2, ...
    'Color', style.dilation, 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top');
text(axB, -4.7, 0.34, 'surrogate band', 'FontName', style.fontName, ...
    'FontSize', 6.2, 'Color', [0.45 0.47 0.50], 'VerticalAlignment', 'bottom');
xlabel(axB, 'Lag from high-gamma peak (s)');
ylabel(axB, 'Pupil change (a.u.)');
title(axB, 'Measure the peri-peak pupil excursion');
axB.XLim = [-5 5];
axB.YLim = [-0.62 0.95];
phg.styleAxes(axB);
phg.addPanelLabel(axB, 'b');

%% (c) The semicontinuous outcome
axC = nexttile(layout);
hold(axC, 'on');
edges = linspace(-3, 3, 41);
continuousPart = [randn(1, 150).*0.55 + 0.75, randn(1, 60).*0.5 - 0.85];
counts = histcounts(continuousPart, edges);
centers = edges(1:end-1) + diff(edges)/2;
bar(axC, centers, counts, 1.0, 'FaceColor', style.lightGray, ...
    'EdgeColor', 'none');

% The zero mass is an order of magnitude taller than the continuous part, so
% plotting both on one linear axis flattens the part that carries the signal.
% The zero bar is clipped at the top of the axis and labelled with its height,
% with a break mark, rather than the axis being rescaled to accommodate it.
zeroMass = 470;
axisTop = max(counts) * 1.55;
bar(axC, 0, axisTop, 0.16, 'FaceColor', style.constriction, ...
    'EdgeColor', 'none');
breakY = axisTop * 0.86;
plot(axC, [-0.14 0.14], [breakY breakY]-axisTop*0.022, '-', ...
    'Color', 'w', 'LineWidth', 1.6);
plot(axC, [-0.14 0.14], [breakY breakY]+axisTop*0.022, '-', ...
    'Color', 'w', 'LineWidth', 1.6);
text(axC, 0.18, axisTop*0.98, ...
    {sprintf('%d exact zeros', zeroMass), 'no suprathreshold run'}, ...
    'FontName', style.fontName, 'FontSize', 6.2, ...
    'Color', style.constriction, 'VerticalAlignment', 'top');
text(axC, 1.35, max(counts)*0.95, {'excursion contacts', 'enter parts 2 and 3'}, ...
    'FontName', style.fontName, 'FontSize', 6.2, 'Color', [0.45 0.47 0.50], ...
    'VerticalAlignment', 'bottom');
xlabel(axC, 'Signed peri-peak area');
ylabel(axC, 'Contacts');
title(axC, 'One outcome, two processes (schematic)');
axC.XLim = [-3 3];
axC.YLim = [0 axisTop];
phg.styleAxes(axC);
phg.addPanelLabel(axC, 'c');

%% (d) The hurdle model
axD = nexttile(layout);
hold(axD, 'on');
axD.XLim = [0 1];
axD.YLim = [0 1];
axis(axD, 'off');
boxes = { ...
    'Part 1 — Prevalence', 'Any excursion at this contact?', ...
        'binomial GLME, all contacts'; ...
    'Part 2 — Direction', 'Given an excursion, dilation or constriction?', ...
        'binomial GLME, excursion contacts (primary)'; ...
    'Part 3 — Magnitude', 'Given an excursion, how large?', ...
        'linear MME on asinh scale, Satterthwaite d.f.'};
boxColor = [0.38 0.41 0.45; style.dilation; style.constriction];
top = 0.96;
height_ = 0.235;
gap = 0.095;
for k = 1:3
    yTop = top - (k-1)*(height_+gap);
    rectangle(axD, 'Position', [0.02 yTop-height_ 0.96 height_], ...
        'Curvature', 0.16, 'FaceColor', [1 1 1], ...
        'EdgeColor', boxColor(k,:), 'LineWidth', 1.0);
    text(axD, 0.06, yTop-0.052, boxes{k,1}, 'FontName', style.fontName, ...
        'FontSize', 7.2, 'FontWeight', 'bold', 'Color', boxColor(k,:));
    text(axD, 0.06, yTop-0.122, boxes{k,2}, 'FontName', style.fontName, ...
        'FontSize', 6.5, 'Color', [0.25 0.27 0.30]);
    text(axD, 0.06, yTop-0.186, boxes{k,3}, 'FontName', style.fontName, ...
        'FontSize', 6.2, 'FontAngle', 'italic', 'Color', [0.45 0.47 0.50]);
    if k < 3
        annotationArrow(axD, 0.50, yTop-height_, yTop-height_-gap+0.012);
    end
end
text(axD, 0.5, 0.035, ...
    'contacts nested in shafts nested in patients; region families FDR-controlled', ...
    'FontName', style.fontName, 'FontSize', 6.2, 'Color', [0.45 0.47 0.50], ...
    'HorizontalAlignment', 'center');
% styleAxes cannot be used on an axis-off panel, so the title is matched to the
% other three by hand; left to itself it renders at the default bold 11 pt and
% dominates the figure.
title(axD, 'Model the two processes separately', ...
    'FontName', style.fontName, 'FontSize', style.titleFontSize, ...
    'FontWeight', 'normal', 'Color', [0.15 0.15 0.15]);
axD.Title.Visible = 'on';
phg.addPanelLabel(axD, 'd');

qc = phg.exportPublicationFigure(fig, cfg, "Fig0_methods_schematic");
close(fig);
end

% -------------------------------------------------------------------------
function smoothed = filterSmooth(signal, window)
kernel = ones(1, window) ./ window;
smoothed = conv(signal, kernel, 'same');
end

function [startIdx, stopIdx] = longestRun(mask)
%LONGESTRUN First and last index of the longest run of true values.
mask = mask(:)';
starts = find(diff([false mask]) == 1);
stops  = find(diff([mask false]) == -1);
[~, best] = max(stops - starts);
startIdx = starts(best);
stopIdx = stops(best);
end

function annotationArrow(ax, x, yFrom, yTo)
plot(ax, [x x], [yFrom yTo], '-', 'Color', [0.60 0.62 0.65], 'LineWidth', 0.8);
plot(ax, x, yTo, 'v', 'MarkerSize', 3.6, 'MarkerFaceColor', [0.60 0.62 0.65], ...
    'MarkerEdgeColor', 'none');
end
