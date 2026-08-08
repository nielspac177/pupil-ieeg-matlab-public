function qc = makeReplicationFigure(cfg, results, options)
%MAKEREPLICATIONFIGURE The replication result, whatever it turned out to be.
%
% Built to be honest about an underpowered test: the primary estimate is shown
% against the discovery interval so that "consistent but wide" is visually
% distinguishable from "replicated", which is the distinction the prespecified
% decision rule turns on and the one a reader is most likely to blur.

arguments
    cfg (1,1) struct
    results (1,1) struct
    options.stem (1,1) string = "FigR2_replication_result"
end

style = phg.figureStyle;
figureHandle = figure('Units', 'inches', 'Position', [1 1 7.2 5.6], ...
    'Color', 'w');
layout = tiledlayout(figureHandle, 2, 2);
phg.setSafeLayout(layout);

%% ------------------------------------------------------- panel a: forest
ax = nexttile(layout);
rows = results.forest;
n = height(rows);
y = n:-1:1;
hold(ax, 'on');

% Discovery interval as a reference band, drawn first so estimates sit on top.
discovery = cfg.replication.discoveryInterval;
patch(ax, [discovery(1) discovery(2) discovery(2) discovery(1)], ...
    [0.4 0.4 n + 0.6 n + 0.6], style.lightGray, 'EdgeColor', 'none', ...
    'FaceAlpha', 0.45);
plot(ax, [1 1], [0.4 n + 0.6], '-', 'Color', [0.35 0.35 0.35], ...
    'LineWidth', 0.75);

for k = 1:n
    % A separated fit is drawn as text, not as an interval: plotting a bar
    % nineteen orders of magnitude wide would imply it means something.
    if ~isfinite(rows.odds_ratio(k)) || rows.separation_suspected(k)
        note = " not estimable";
        if rows.separation_suspected(k)
            note = sprintf(' not estimable (interval spans %.0f decades)', ...
                rows.ci_width_decades(k));
        end
        text(ax, 1, y(k), note, 'FontName', style.fontName, ...
            'FontSize', style.tickFontSize, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'Color', [0.4 0.4 0.4]);
        continue
    end
    isPrimary = rows.is_primary(k);
    colour = style.gray;
    if isPrimary
        colour = style.constriction;
    end
    plot(ax, [rows.ci_low(k) rows.ci_high(k)], [y(k) y(k)], '-', ...
        'Color', colour, 'LineWidth', style.lineWidth);
    plot(ax, rows.odds_ratio(k), y(k), 'o', 'MarkerSize', 4.5, ...
        'MarkerFaceColor', colour, 'MarkerEdgeColor', 'none');
end
hold(ax, 'off');

% The x range is set from the estimable intervals and the discovery band only.
% Letting a separated fit into the limits stretched the axis across twenty
% decades and squeezed the primary estimate into a few pixels.
estimable = isfinite(rows.odds_ratio) & ~rows.separation_suspected;
if any(estimable)
    lower = min([rows.ci_low(estimable); discovery(1)]);
    upper = max([rows.ci_high(estimable); discovery(2); 1]);
else
    lower = discovery(1);
    upper = 1;
end
xlim(ax, [lower / 2, upper * 2]);
set(ax, 'XScale', 'log', 'YTick', y(end:-1:1), ...
    'YTickLabel', flipud(rows.label));
ylim(ax, [0.4, n + 0.6]);
xlabel(ax, 'Odds of dilation-linkage, hippocampal vs extrahippocampal');
title(ax, 'H1 estimates vs discovery interval (shaded)');
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'a');

%% -------------------------------------------- panel b: per-subject rates
ax = nexttile(layout);
paired = results.perSubject;
hold(ax, 'on');
for k = 1:height(paired)
    if ~isfinite(paired.hippocampal_rate(k))
        continue
    end
    plot(ax, [1 2], [paired.extrahippocampal_rate(k), ...
        paired.hippocampal_rate(k)], '-o', 'Color', style.gray, ...
        'MarkerSize', 4, 'MarkerFaceColor', 'w', 'LineWidth', 0.9);
end
hold(ax, 'off');
xlim(ax, [0.7 2.3]);
ylim(ax, [-0.05 1.05]);
set(ax, 'XTick', [1 2], 'XTickLabel', {'Extrahippocampal', 'Hippocampal'});
ylabel(ax, 'Proportion dilation-linked');
title(ax, sprintf('Per participant (%d with hippocampal coverage)', ...
    sum(isfinite(paired.hippocampal_rate))));
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'b');

%% ------------------------------------------ panel c: cohort comparison
ax = nexttile(layout);
comparison = results.cohortRates;
bars = bar(ax, comparison.values, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
    'BarWidth', 0.75);
bars.CData(1, :) = style.dilation;
bars.CData(2, :) = style.constriction;
set(ax, 'XTickLabel', comparison.labels);
ylabel(ax, 'Proportion dilation-linked');
ylim(ax, [0 1.05]);
for k = 1:numel(comparison.values)
    text(ax, k, comparison.values(k) + 0.03, ...
        sprintf('%.0f%%\n(n=%d)', 100 * comparison.values(k), ...
        comparison.counts(k)), 'HorizontalAlignment', 'center', ...
        'FontName', style.fontName, 'FontSize', style.tickFontSize);
end
% Named as unadjusted on purpose: these raw proportions run opposite to the
% mixed-model estimate in panel a, because the model separates within- from
% between-participant variation and these do not.
title(ax, 'Unadjusted proportions (cf. model in a)');
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'c');

%% ------------------------------------------------ panel d: ripple test
ax = nexttile(layout);
if isfield(results, 'ripple') && ~isempty(results.ripple) && ...
        isfinite(results.ripple.observed)
    histogram(ax, results.ripple.null, 24, 'FaceColor', style.lightGray, ...
        'EdgeColor', 'none', 'Normalization', 'probability');
    yLimits = ylim(ax);
    hold(ax, 'on');
    plot(ax, [results.ripple.observed results.ripple.observed], yLimits, ...
        '-', 'Color', style.dilation, 'LineWidth', 1.6);
    hold(ax, 'off');
    ylim(ax, yLimits);
    xlabel(ax, 'Fraction of ripples during pupil constriction');
    ylabel(ax, 'Circular-shift surrogates');
    title(ax, sprintf('H2: %d ripples, observed %.1f%%, P = %.3g', ...
        results.ripple.nRipples, 100 * results.ripple.observed, ...
        results.ripple.pValue));
else
    text(ax, 0.5, 0.5, {'H2 not estimable', ...
        'no accepted ripples in eligible epochs'}, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'FontName', style.fontName, ...
        'FontSize', style.labelFontSize, 'Color', [0.4 0.4 0.4]);
    set(ax, 'XTick', [], 'YTick', []);
    title(ax, 'H2: hippocampal ripples and pupil state');
end
phg.styleAxes(ax);
phg.addPanelLabel(ax, 'd');

qc = phg.exportPublicationFigure(figureHandle, cfg, options.stem);
close(figureHandle);
end
