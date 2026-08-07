function qc = auditFigureLayout(figureHandle, imagePath, stem, options)
%AUDITFIGURELAYOUT Check a figure for clipping, dead margins, and text overlap.
%
%   The previous audit only asked whether rendered content came within four
%   pixels of the export boundary. That test passes trivially on a figure that
%   is one third blank canvas, and it cannot see label collisions at all, so it
%   certified figures that had exactly the problems it was meant to catch.
%
%   This audit reports three independent conditions:
%     clipping      content touching or nearly touching an export boundary
%     dead margin   a blank border wider than MaxMarginFraction of that
%                   dimension, measured from the rendered raster
%     overlap       pairwise intersection of text bounding boxes, measured in
%                   figure pixel coordinates before the figure is closed
%
%   A figure passes only if all three are clear.

arguments
    figureHandle (1,1) matlab.ui.Figure
    imagePath (1,1) string
    stem (1,1) string
    options.MaxMarginCentimetres (1,1) double = 1.2
    options.MinMarginPixels (1,1) double = 4
    options.MaxOverlapFraction (1,1) double = 0.08
end

%% Raster margins
imageData = imread(imagePath);
if ismatrix(imageData)
    nonwhite = imageData < 250;
else
    nonwhite = any(imageData(:,:,1:min(3,size(imageData,3))) < 250, 3);
end
rows = find(any(nonwhite, 2));
columns = find(any(nonwhite, 1));
height_ = size(nonwhite, 1);
width_ = size(nonwhite, 2);

% A dead margin is judged in printed centimetres rather than as a fraction of
% the canvas: what matters to a reader is how much blank paper surrounds the
% figure at its published size, not how that compares to the figure's own area.
originalUnits = figureHandle.Units;
figureHandle.Units = 'centimeters';
figureWidthCm = figureHandle.Position(3);
figureHandle.Units = originalUnits;
pixelsPerCm = width_ / figureWidthCm;

if isempty(rows) || isempty(columns)
    margins = [NaN NaN NaN NaN];
    boundaryClear = false;
    marginsTight = false;
    whitespaceFraction = 1;
else
    margins = [rows(1)-1, height_-rows(end), columns(1)-1, width_-columns(end)];
    boundaryClear = all(margins >= options.MinMarginPixels);
    marginsTight = all(margins ./ pixelsPerCm <= options.MaxMarginCentimetres);
    contentArea = (rows(end)-rows(1)+1) * (columns(end)-columns(1)+1);
    whitespaceFraction = 1 - contentArea / (height_ * width_);
end

%% Text overlap in figure pixel coordinates
[overlapCount, worstOverlap, worstPair] = localTextOverlap(figureHandle, ...
    options.MaxOverlapFraction);

layoutOk = boundaryClear && marginsTight && overlapCount == 0;

qc = table(stem, width_, height_, margins(1), margins(2), margins(3), ...
    margins(4), max(margins) / pixelsPerCm, whitespaceFraction, ...
    boundaryClear, marginsTight, overlapCount, worstOverlap, ...
    string(worstPair), layoutOk, 'VariableNames', {'figure', ...
    'width_pixels', 'height_pixels', 'top_margin_pixels', ...
    'bottom_margin_pixels', 'left_margin_pixels', 'right_margin_pixels', ...
    'largest_margin_cm', 'whitespace_fraction', 'boundary_clear', ...
    'margins_tight', 'text_overlap_count', 'worst_overlap_fraction', ...
    'worst_overlap_pair', 'layout_ok'});
end

% -------------------------------------------------------------------------
function [overlapCount, worstOverlap, worstPair] = localTextOverlap(fig, tolerance)
%LOCALTEXTOVERLAP Pairwise intersection of visible text boxes, in figure pixels.

boxes = zeros(0, 4);
labels = strings(0, 1);

% Annotation textboxes live in a hidden figure-level axes, not in any data
% axes, so the per-axes sweep below cannot see them. They are collected first
% in figure-normalized units and converted to pixels.
figurePixels = getpixelposition(fig);
annotations = findall(fig, 'Type', 'textboxshape');
for k = 1:numel(annotations)
    a = annotations(k);
    if strcmp(a.Visible, 'off') || strlength(strtrim(join(string(a.String), ""))) == 0
        continue
    end
    position = a.Position;
    boxes(end+1, :) = [position(1)*figurePixels(3), ...
        position(2)*figurePixels(4), position(3)*figurePixels(3), ...
        position(4)*figurePixels(4)]; %#ok<AGROW>
    text_ = string(a.String);
    labels(end+1, 1) = strtrim(text_(1)); %#ok<AGROW>
end

allAxes = findall(fig, 'Type', 'axes');
for a = 1:numel(allAxes)
    ax = allAxes(a);
    if strcmp(ax.Visible, 'off') && isempty(findall(ax, 'Type', 'text'))
        continue
    end
    axPixels = getpixelposition(ax, true);
    % findall searches hidden handles, so it already returns Title, XLabel and
    % YLabel. Appending them again would make every one of them overlap itself.
    candidates = findall(ax, 'Type', 'text');
    for k = 1:numel(candidates)
        t = candidates(k);
        if ~isgraphics(t) || strcmp(t.Visible, 'off') || isempty(char(join(string(t.String), "")))
            continue
        end
        originalUnits = t.Units;
        try
            t.Units = 'pixels';
            extent = t.Extent;
            t.Units = originalUnits;
        catch
            t.Units = originalUnits;
            continue
        end
        if any(~isfinite(extent)) || extent(3) <= 0 || extent(4) <= 0
            continue
        end
        boxes(end+1, :) = [axPixels(1)+extent(1), axPixels(2)+extent(2), ...
            extent(3), extent(4)]; %#ok<AGROW>
        text_ = string(t.String);
        labels(end+1, 1) = strtrim(text_(1)); %#ok<AGROW>
    end
end

overlapCount = 0;
worstOverlap = 0;
worstPair = "none";
for i = 1:size(boxes, 1)
    for j = i+1:size(boxes, 1)
        overlapWidth = min(boxes(i,1)+boxes(i,3), boxes(j,1)+boxes(j,3)) - ...
            max(boxes(i,1), boxes(j,1));
        overlapHeight = min(boxes(i,2)+boxes(i,4), boxes(j,2)+boxes(j,4)) - ...
            max(boxes(i,2), boxes(j,2));
        if overlapWidth <= 0 || overlapHeight <= 0
            continue
        end
        intersection = overlapWidth * overlapHeight;
        smaller = min(boxes(i,3)*boxes(i,4), boxes(j,3)*boxes(j,4));
        fraction = intersection / smaller;
        if fraction > tolerance
            overlapCount = overlapCount + 1;
            if fraction > worstOverlap
                worstOverlap = fraction;
                worstPair = labels(i) + " | " + labels(j);
            end
        end
    end
end
end
