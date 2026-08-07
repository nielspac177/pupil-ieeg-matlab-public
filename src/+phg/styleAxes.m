function styleAxes(ax)
%STYLEAXES Apply consistent manuscript-scale axes styling.

style = phg.figureStyle;
set(ax, 'FontName', style.fontName, 'FontSize', style.tickFontSize, ...
    'LineWidth', style.axisLineWidth, 'TickDir', 'out', ...
    'Box', 'off', 'Color', 'w', 'Layer', 'top');
ax.XLabel.FontSize = style.labelFontSize;
ax.YLabel.FontSize = style.labelFontSize;
ax.Title.FontSize = style.titleFontSize;
ax.Title.FontWeight = 'normal';
end
