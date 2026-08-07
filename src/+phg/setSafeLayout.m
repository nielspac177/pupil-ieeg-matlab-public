function setSafeLayout(layout, options)
%SETSAFELAYOUT Give a tiled layout tight, content-driven margins.
%   The previous implementation pinned every layout to a fixed OuterPosition of
%   [0.06 0.10 0.88 0.78], reserving 12% of the height and 6% of the width as
%   blank canvas on every figure regardless of its content. At 600 dpi that is
%   roughly one blank inch per side.
%
%   MATLAB's 'compact' and 'tight' modes size margins from the decorations that
%   are actually present, so the layout is allowed to fill the figure and only
%   the strip needed for panel labels is reserved.

arguments
    layout (1,1) matlab.graphics.layout.TiledChartLayout
    options.TileSpacing (1,1) string = "compact"
    options.Padding (1,1) string = "compact"
    options.PanelLabelRoom (1,1) logical = true
end

layout.TileSpacing = char(options.TileSpacing);
layout.Padding = char(options.Padding);
layout.Units = 'normalized';

if options.PanelLabelRoom
    % phg.addPanelLabel writes at x = -0.04, y = 1.03 in axes-normalized units,
    % so a narrow strip on the left and top is reserved for it.
    layout.OuterPosition = [0.015 0.005 0.980 0.965];
else
    layout.OuterPosition = [0 0 1 1];
end
end
