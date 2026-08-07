function figureHandle = makeElectrodeLocalization(T, cfg, options)
%MAKEELECTRODELOCALIZATION MNI electrode coverage figures.
%
%   Two renderers are available.
%
%   "legui3d" (default) reproduces the rendering LeGUI itself uses: the
%   template surface produced by LeG_genSurfaces, drawn as a semi-transparent
%   grey patch with gouraud lighting and LeGUI's own reflectance settings
%   (FaceColor [0.5 0.5 0.5], ambient 0.3, diffuse 0.8, specular 0.1) lit by a
%   headlight, with contacts as three-dimensional markers inside it.
%
%   "glass" projects the same surface to a two-dimensional silhouette. It
%   resolves overlapping depth contacts better but is not what LeGUI draws.
%
%   The earlier 3-D implementation failed on layout, not on rendering: a 2x2
%   grid with a full-height colour bar left the brains occupying a small
%   fraction of the canvas. Here the columns are sized from each view's
%   projected aspect ratio and the colour key sits in a short south strip.
%
%   The overlay set is contacts with a suprathreshold pupil excursion
%   (RespAreaNet ~= 0), which is threshold-free, rather than the historical
%   RespSig > 0.10 selection.

arguments
    T table
    cfg (1,1) struct
    options.mode (1,1) string {mustBeMember(options.mode, ...
        ["effect", "region", "polarity"])} = "effect"
    options.renderer (1,1) string {mustBeMember(options.renderer, ...
        ["legui3d", "glass"])} = "legui3d"
end

style = phg.figureStyle;
xyz = T.XYZMNI;
valid = all(isfinite(xyz), 2);
overlay = T.RespAreaNet ~= 0 & valid;
region = phg.cleanRegionLabels(T.NMM);
signedEffect = T.RespAreaNet;

surface = phg.buildTransparentBrainSurface(cfg);

% Screen bases for each panel, expressed in MNI axes. Signs follow the
% orientation a reader expects: anterior away from the midline on the lateral
% views, anterior upward on the axial view. The camera fields are the
% equivalent azimuth/elevation for the three-dimensional renderer.
views = struct( ...
    'name', {'Left lateral', 'Anterior', 'Right lateral', 'Superior'}, ...
    'u', {[0 -1 0], [-1 0 0], [0 1 0], [1 0 0]}, ...
    'v', {[0 0 1], [0 0 1], [0 0 1], [0 1 0]}, ...
    'camera', {[-90 0], [180 0], [90 0], [0 90]}, ...
    'hemisphere', {-1, 0, 1, 0}, ...
    'leftLabel', {'', 'R', '', 'L'}, ...
    'rightLabel', {'', 'L', '', 'R'});

[overlayColors, legendLabels, legendColors, colorValues, effectLimit, stem] = ...
    configureMode(options.mode);
if options.renderer == "glass"
    stem = replace(stem, "_legui", "_glass");
end

% Each view projects the MNI box to a different aspect ratio: the lateral views
% are wide and short (anterior-posterior by superior-inferior), the axial view
% is narrow and long. Equal tiles would leave the wide views width-limited and
% the axial view flanked by blank canvas, so column width is allocated in
% proportion to each view's projected aspect using a finer underlying grid.
% This applies to both renderers because the 3-D renderer sets the data aspect
% with daspect rather than "axis equal", and so keeps the rectangle the layout
% assigned it instead of shrinking to a square.
columnSpans = [6 5 6 4];
totalColumns = sum(columnSpans);
columnStart = cumsum([1 columnSpans(1:end-1)]);

figureHandle = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 18.3 5.4]);
layout = tiledlayout(figureHandle, 1, totalColumns, ...
    'TileSpacing', 'tight', 'Padding', 'compact');
phg.setSafeLayout(layout, 'TileSpacing', "tight");

panelLabels = ["a", "b", "c", "d"];
axesHandles = gobjects(4, 1);
titleHandles = gobjects(4, 1);
labelHandles = gobjects(4, 1);
for k = 1:4
    axesHandles(k) = nexttile(layout, columnStart(k), [1 columnSpans(k)]);
    if options.renderer == "legui3d"
        plotLeGUIView(axesHandles(k), views(k), panelLabels(k));
    else
        plotGlassView(axesHandles(k), views(k), panelLabels(k));
    end
    titleHandles(k) = axesHandles(k).Title;
    labelHandles(k) = findobj(axesHandles(k), 'Tag', 'phgPanelLabel');
end

if options.mode == "effect"
    colormap(figureHandle, phg.divergingMap(257));
    for k = 1:4
        clim(axesHandles(k), [-effectLimit effectLimit]);
    end
    % Positioned in figure coordinates rather than given its own layout tile.
    % A tile-managed colour bar takes the full width of the strip and a fixed
    % share of its height, which on a short four-panel figure renders a bar
    % wider and thicker than the brains it annotates.
    colorbarHandle = colorbar(axesHandles(4));
    drawnow;
    colorbarHandle.Units = 'normalized';
    colorbarHandle.Location = 'south';
    colorbarHandle.Position = [0.335 0.085 0.33 0.045];
    colorbarHandle.AxisLocation = 'out';
    colorbarHandle.FontName = style.fontName;
    colorbarHandle.FontSize = style.tickFontSize;
    % The centre tick is drawn but left unlabelled: the caption sits under the
    % middle of the bar, and a "+0" there would collide with it. The zero point
    % of a diverging map is legible from the white midpoint anyway.
    colorbarHandle.Ticks = linspace(-effectLimit, effectLimit, 3);
    endLabels = cellstr(compose('%+.0f', colorbarHandle.Ticks));
    endLabels{2} = '';
    colorbarHandle.TickLabels = endLabels;
    % The colour bar's own Label is not laid out reliably once Position is set
    % by hand, so the caption is drawn as a figure-level annotation instead.
    annotation(figureHandle, 'textbox', [0.235 0.005 0.53 0.06], ...
        'String', 'Signed pupil-coupling area (a.u.)', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', 'FontName', style.fontName, ...
        'FontSize', style.labelFontSize);
else
    legendHandles = gobjects(numel(legendLabels), 1);
    for k = 1:numel(legendLabels)
        legendHandles(k) = scatter(axesHandles(1), NaN, NaN, 26, ...
            legendColors(k,:), 'filled', 'MarkerEdgeColor', 'white');
    end
    legendHandle = legend(axesHandles(1), legendHandles, ...
        cellstr(legendLabels), 'Orientation', 'horizontal', ...
        'NumColumns', numel(legendLabels), ...
        'FontName', style.fontName, 'FontSize', 6.3, 'Box', 'off');
    legendHandle.Layout.Tile = 'south';
end

phg.exportPublicationFigure(figureHandle, cfg, stem, ...
    'pdfContentType', "vector");
fprintf("[PHG] Wrote %s localization figure (%s renderer).\n", ...
    options.mode, options.renderer);

    function plotLeGUIView(ax, viewSpec, panelLabel)
        %PLOTLEGUIVIEW Render one view the way LeGUI renders its brain surface.
        hold(ax, 'on');

        if ~isempty(surface.vertices)
            patch(ax, 'Faces', surface.faces, 'Vertices', surface.vertices, ...
                'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', cfg.brainAlpha, ...
                'AlphaDataMapping', 'none', 'EdgeColor', 'none', ...
                'FaceLighting', 'gouraud', 'AmbientStrength', 0.3, ...
                'DiffuseStrength', 0.8, 'SpecularStrength', 0.1, ...
                'SpecularExponent', 10);
        end

        % Both hemispheres are shown in every view. A lateral view of a
        % bilateral implant would otherwise hide half the contacts behind the
        % midline, and the transparent surface already conveys depth.
        background = valid & ~overlay;
        foreground = valid & overlay;

        scatter3(ax, xyz(background,1), xyz(background,2), xyz(background,3), ...
            7, [0.66 0.68 0.71], 'filled', 'MarkerFaceAlpha', 0.30, ...
            'MarkerEdgeAlpha', 0);
        if options.mode == "effect"
            scatter3(ax, xyz(foreground,1), xyz(foreground,2), ...
                xyz(foreground,3), 26, colorValues(foreground), 'filled', ...
                'MarkerEdgeColor', [1 1 1], 'LineWidth', 0.25, ...
                'MarkerFaceAlpha', 0.95);
        else
            scatter3(ax, xyz(foreground,1), xyz(foreground,2), ...
                xyz(foreground,3), 26, overlayColors(foreground,:), ...
                'filled', 'MarkerEdgeColor', [1 1 1], 'LineWidth', 0.25, ...
                'MarkerFaceAlpha', 0.95);
        end

        points = xyz(valid,:);
        if ~isempty(surface.vertices)
            points = [points; surface.vertices];
        end
        bounds = [min(points, [], 1); max(points, [], 1)];
        padding = max(diff(bounds, 1, 1) .* 0.02, 1);
        xlim(ax, bounds(:,1)' + [-padding(1) padding(1)]);
        ylim(ax, bounds(:,2)' + [-padding(2) padding(2)]);
        zlim(ax, bounds(:,3)' + [-padding(3) padding(3)]);

        % daspect rather than "axis equal": both give isotropic data units, but
        % "axis equal" also resizes each axes rectangle to the content, which
        % leaves the four panels at different heights and their titles and
        % panel labels on different lines. Setting the data aspect directly
        % keeps every axes at the rectangle the tiled layout assigned it.
        daspect(ax, [1 1 1]);
        axis(ax, 'vis3d', 'off');
        view(ax, viewSpec.camera);
        camproj(ax, 'orthographic');
        camlight(ax, 'headlight');
        % With "axis vis3d" MATLAB fits the bounding sphere of the data box,
        % not its projection, so a brain occupies only about 1/1.5 of the axes
        % by default. The zoom recovers the difference between the sphere
        % diameter and the widest projected extent.
        boxSpan = diff(bounds, 1, 1);
        sphereDiameter = norm(boxSpan);
        projectedSpan = max(abs(boxSpan * viewSpec.u(:)), ...
            abs(boxSpan * viewSpec.v(:)));
        % The zoom factor is bounded above by content clipping, not by the axes
        % box: camzoom moves the camera without changing the data limits, so
        % MATLAB does not clip the projection to the tile rectangle. Some
        % contacts sit outside the template surface, and above about 0.80 they
        % reach the export boundary. 0.78 is the largest value at which all
        % three modes export with a clear boundary; the residual side margin of
        % roughly 1.5 cm is the cost of that headroom.
        camzoom(ax, 0.78 * sphereDiameter / projectedSpan);

        if ~isempty(viewSpec.leftLabel)
            text(ax, 0.05, 0.13, viewSpec.leftLabel, 'Units', 'normalized', ...
                'FontName', style.fontName, 'FontSize', 7, ...
                'Color', [0.42 0.45 0.48], 'HorizontalAlignment', 'left');
            text(ax, 0.95, 0.13, viewSpec.rightLabel, 'Units', 'normalized', ...
                'FontName', style.fontName, 'FontSize', 7, ...
                'Color', [0.42 0.45 0.48], 'HorizontalAlignment', 'right');
        end
        title(ax, viewSpec.name, 'FontName', style.fontName, ...
            'FontSize', style.titleFontSize, 'FontWeight', 'normal');
        phg.addPanelLabel(ax, panelLabel);
        set(ax, 'Color', 'w');
    end

    function plotGlassView(ax, viewSpec, panelLabel)
        %PLOTGLASSVIEW Render one view as a flat projected silhouette.
        hold(ax, 'on');

        projection = struct('u', viewSpec.u, 'v', viewSpec.v);
        if viewSpec.hemisphere ~= 0
            hemisphere = viewSpec.hemisphere;
            projection.mask = @(V) sign(hemisphere) * V(:,1) >= -1;
        end
        outline = phg.glassBrainOutline(surface, projection);
        if ~isempty(outline.u)
            fill(ax, outline.u, outline.v, [0.94 0.945 0.95], ...
                'EdgeColor', [0.62 0.65 0.68], 'LineWidth', 0.7);
        end
        if ~isempty(outline.midlineU)
            plot(ax, outline.midlineU, outline.midlineV, '-', ...
                'Color', [0.80 0.83 0.86], 'LineWidth', 0.55);
        end

        inView = valid;
        if viewSpec.hemisphere < 0
            inView = inView & xyz(:,1) <= 1;
        elseif viewSpec.hemisphere > 0
            inView = inView & xyz(:,1) >= -1;
        end
        background = inView & ~overlay;
        foreground = inView & overlay;

        u = xyz * viewSpec.u(:);
        v = xyz * viewSpec.v(:);

        scatter(ax, u(background), v(background), 5, [0.66 0.68 0.71], ...
            'filled', 'MarkerFaceAlpha', 0.40, 'MarkerEdgeAlpha', 0);
        if options.mode == "effect"
            scatter(ax, u(foreground), v(foreground), 22, ...
                colorValues(foreground), 'filled', ...
                'MarkerEdgeColor', [1 1 1], 'LineWidth', 0.25, ...
                'MarkerFaceAlpha', 0.92);
        else
            scatter(ax, u(foreground), v(foreground), 22, ...
                overlayColors(foreground,:), 'filled', ...
                'MarkerEdgeColor', [1 1 1], 'LineWidth', 0.25, ...
                'MarkerFaceAlpha', 0.92);
        end

        axis(ax, 'equal');
        axis(ax, 'off');
        limits = [min([outline.u; u(inView)]), max([outline.u; u(inView)]), ...
            min([outline.v; v(inView)]), max([outline.v; v(inView)])];
        padding = 0.04 * max(limits(2)-limits(1), limits(4)-limits(3));
        xlim(ax, [limits(1)-padding limits(2)+padding]);
        ylim(ax, [limits(3)-padding limits(4)+padding*2.2]);

        if ~isempty(viewSpec.leftLabel)
            text(ax, 0.04, 0.90, viewSpec.leftLabel, 'Units', 'normalized', ...
                'FontName', style.fontName, 'FontSize', 7, ...
                'Color', [0.42 0.45 0.48], 'HorizontalAlignment', 'left');
            text(ax, 0.96, 0.90, viewSpec.rightLabel, 'Units', 'normalized', ...
                'FontName', style.fontName, 'FontSize', 7, ...
                'Color', [0.42 0.45 0.48], 'HorizontalAlignment', 'right');
        end
        title(ax, viewSpec.name, 'FontName', style.fontName, ...
            'FontSize', style.titleFontSize, 'FontWeight', 'normal');
        phg.addPanelLabel(ax, panelLabel);
        set(ax, 'Color', 'w');
    end

    function [colors, labels, keyColors, values, limit, outputStem] = ...
            configureMode(mode)
        colors = repmat(style.gray, height(T), 1);
        values = signedEffect;
        limit = 1;
        switch mode
            case "effect"
                finiteEffect = signedEffect(overlay & isfinite(signedEffect));
                limit = prctile(abs(finiteEffect), 95);
                if isempty(limit) || ~isfinite(limit) || limit <= 0
                    limit = 1;
                end
                labels = strings(0,1);
                keyColors = zeros(0,3);
                outputStem = "Fig7_electrode_effect_legui";
            case "region"
                overlayRegions = region(overlay);
                [names, ~, index] = unique(overlayRegions);
                counts = accumarray(index, 1);
                [~, order] = sort(counts, 'descend');
                % Four labelled regions, not five. The fifth-ranked region
                % (MFG) is assigned a dark grey by phg.regionColor, which on a
                % grey brain is not separable from either the pooled "Other"
                % grey or the pale "all contacts" grey. Keeping only regions
                % with distinct hues leaves exactly two greys in the legend,
                % differing clearly in lightness.
                topNames = names(order(1:min(4, numel(order))));
                colors = repmat([0.42 0.44 0.47], height(T), 1);
                for r = 1:numel(topNames)
                    mask = region == topNames(r);
                    colors(mask,:) = repmat(phg.regionColor(topNames(r)), ...
                        sum(mask), 1);
                end
                labels = string(sprintf("All contacts (n=%d)", sum(valid)));
                keyColors = [0.66 0.68 0.71];
                for r = 1:numel(topNames)
                    labels(end+1,1) = string(sprintf("%s (n=%d)", ...
                        char(phg.shortRegionLabel(topNames(r))), ...
                        sum(overlay & region == topNames(r)))); %#ok<AGROW>
                    keyColors(end+1,:) = phg.regionColor(topNames(r)); %#ok<AGROW>
                end
                otherCount = sum(overlay & ~ismember(region, topNames));
                labels(end+1,1) = string(sprintf("Other (n=%d)", otherCount));
                keyColors(end+1,:) = [0.42 0.44 0.47];
                outputStem = "Fig8_electrode_regions_legui";
            case "polarity"
                dilation = overlay & signedEffect > 0;
                constriction = overlay & signedEffect < 0;
                colors(dilation,:) = repmat(style.dilation, sum(dilation), 1);
                colors(constriction,:) = repmat(style.constriction, ...
                    sum(constriction), 1);
                labels = [string(sprintf("All contacts (n=%d)", sum(valid))); ...
                    string(sprintf("Dilation-linked (n=%d)", sum(dilation))); ...
                    string(sprintf("Constriction-linked (n=%d)", ...
                    sum(constriction)))];
                keyColors = [0.66 0.68 0.71; style.dilation; style.constriction];
                outputStem = "Fig9_electrode_polarity_legui";
        end
    end
end
