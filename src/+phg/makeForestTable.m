function qc = makeForestTable(cfg, analysis)
%MAKEFORESTTABLE Region results as a JAMA-style table with an inline forest.
%
%   A table and a forest plot answer different questions — the table gives the
%   reader an exact number to quote, the forest gives them the comparison at a
%   glance — and journals increasingly print both in one object rather than
%   making the reader hold a table and a figure in mind at once.
%
%   Two axes share a row grid: the left carries the text columns, the right the
%   estimate and its interval on a logarithmic scale with a reference line at
%   unity. Rows are ordered by effect size, and rows surviving FDR correction
%   are drawn filled so significance is legible without reading the q column.

arguments
    cfg (1,1) struct
    analysis (1,1) struct
end

style = phg.figureStyle;

prevalence = analysis.prevalenceTable;
prevalence = prevalence(prevalence.in_fdr_family, :);
prevalence.region = extractAfter(prevalence.term, "Region_");
prevalence = sortrows(prevalence, 'odds_ratio', 'descend');

coverage = readtable(fullfile(cfg.tableDir, 'region_coverage_summary.csv'), ...
    'TextType', 'string');
[~, position] = ismember(prevalence.region, coverage.region);

nRows = height(prevalence);
fig = figure('Color', 'w', 'Units', 'centimeters', 'Visible', 'off', ...
    'Position', [2 2 18.0 1.30 + 0.52*nRows]);
layout = tiledlayout(fig, 1, 12, 'TileSpacing', 'none', 'Padding', 'compact');
layout.OuterPosition = [0.010 0.010 0.980 0.965];

%% Text columns
axText = nexttile(layout, 1, [1 8]);
hold(axText, 'on');
axis(axText, 'off');
xlim(axText, [0 1]);
ylim(axText, [0.4 nRows + 1.25]);
set(axText, 'YDir', 'reverse');

% Right edges, not left starts: every column but the label is right-aligned,
% and the intervals are long enough that left-anchored columns collide.
columnEdge = [0.00 0.34 0.58 0.88 1.00];
headers = {'Region', 'Contacts', 'Excursion', 'OR [95% CI]', 'q'};
for c = 1:numel(columnEdge)
    if c == 1
        alignment = 'left';
    else
        alignment = 'right';
    end
    text(axText, columnEdge(c), 0.62, headers{c}, ...
        'FontName', style.fontName, 'FontSize', 7.2, 'FontWeight', 'bold', ...
        'HorizontalAlignment', alignment, 'VerticalAlignment', 'bottom');
end

plot(axText, [0 1], [0.78 0.78], '-', 'Color', [0.25 0.25 0.25], ...
    'LineWidth', 0.9, 'Clipping', 'off');

for k = 1:nRows
    row = prevalence(k, :);
    cov = coverage(position(k), :);
    weight = 'normal';
    if row.fdr_q_value < 0.05
        weight = 'bold';
    end
    text(axText, columnEdge(1), k, phg.shortRegionLabel(row.region), ...
        'FontName', style.fontName, 'FontSize', 7.0, 'FontWeight', weight, ...
        'VerticalAlignment', 'middle');
    text(axText, columnEdge(2), k, string(cov.n_contacts), ...
        'FontName', style.fontName, 'FontSize', 7.0, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    text(axText, columnEdge(3), k, ...
        sprintf('%d (%.0f%%)', cov.n_excursion_contacts, ...
        100*cov.excursion_proportion), 'FontName', style.fontName, ...
        'FontSize', 7.0, 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle');
    text(axText, columnEdge(4), k, ...
        sprintf('%.2f [%.2f, %.2f]', row.odds_ratio, ...
        row.odds_ratio_ci95_low, row.odds_ratio_ci95_high), ...
        'FontName', style.fontName, 'FontSize', 7.0, 'FontWeight', weight, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    if row.fdr_q_value < 0.001
        qLabel = sprintf('%.1e', row.fdr_q_value);
    else
        qLabel = sprintf('%.3f', row.fdr_q_value);
    end
    text(axText, columnEdge(5), k, qLabel, 'FontName', style.fontName, ...
        'FontSize', 7.0, 'FontWeight', weight, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
end

%% Forest
axForest = nexttile(layout, 9, [1 4]);
hold(axForest, 'on');
ylim(axForest, [0.4 nRows + 1.25]);
set(axForest, 'YDir', 'reverse', 'YTick', []);
set(axForest, 'XScale', 'log');
xlim(axForest, [0.05 40]);

xline(axForest, 1, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
for k = 1:nRows
    row = prevalence(k, :);
    colour = phg.regionColor(row.region);
    significant = row.fdr_q_value < 0.05;
    plot(axForest, [row.odds_ratio_ci95_low row.odds_ratio_ci95_high], ...
        [k k], '-', 'Color', colour, 'LineWidth', 1.2);
    % Interval end caps, as in a printed forest plot.
    for edge = [row.odds_ratio_ci95_low row.odds_ratio_ci95_high]
        plot(axForest, [edge edge], k + [-0.16 0.16], '-', ...
            'Color', colour, 'LineWidth', 1.0);
    end
    if significant
        faceColour = colour;
    else
        faceColour = 'white';
    end
    scatter(axForest, row.odds_ratio, k, 30, 'o', 'filled', ...
        'MarkerFaceColor', faceColour, 'MarkerEdgeColor', colour, ...
        'LineWidth', 0.9);
end

axForest.XTick = [0.1 0.5 1 5 20];
axForest.XTickLabel = {'0.1', '0.5', '1', '5', '20'};
xlabel(axForest, 'Odds of an excursion vs rest of implant');
phg.styleAxes(axForest);
axForest.YAxis.Visible = 'off';
axForest.XAxisLocation = 'bottom';

text(axForest, 0.45, 0.62, 'less likely', 'FontName', style.fontName, ...
    'FontSize', 6.2, 'Color', [0.45 0.45 0.45], ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
text(axForest, 2.2, 0.62, 'more likely', 'FontName', style.fontName, ...
    'FontSize', 6.2, 'Color', [0.45 0.45 0.45], ...
    'VerticalAlignment', 'bottom');

qc = phg.exportPublicationFigure(fig, cfg, "Fig11_prevalence_forest_table");
close(fig);
end
