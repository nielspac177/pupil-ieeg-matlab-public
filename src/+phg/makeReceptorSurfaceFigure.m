function qc = makeReceptorSurfaceFigure(channel, cfg)
%MAKERECEPTORSURFACEFIGURE Cortical surface, subcortical meshes, and leads.
%
%   A Lead-DBS-style render: a semi-transparent cortical surface carrying the
%   normative norepinephrine-transporter map, solid hippocampus and amygdala
%   meshes built from the Neuromorphometrics atlas that LeGUI ships, and the
%   implanted electrodes drawn as leads rather than as a scatter of points —
%   contacts on one shaft are joined in order, which is what they physically
%   are.
%
%   Panel a shows the receptor map alone, panel b the same view with the
%   contacts coloured by coupling direction, so a reader can see the gradient
%   and the sampling separately before seeing them together.
%
%   Requires the PET map in <workDir>/cache/pet and the LeGUI atlas at
%   cfg.leguiRoot. Returns empty and says why if either is missing.

arguments
    channel table
    cfg (1,1) struct
end

qc = table.empty;
style = phg.figureStyle;

petFile = fullfile(cfg.workDir, 'cache', 'pet', 'NAT_MRB_hc77_ding.nii');
if ~isfield(cfg, 'leGUIRoot') || cfg.leGUIRoot == ""
    fprintf('[PHG] Receptor surface figure skipped: cfg.leGUIRoot not set.\n');
    return
end
atlasFile = fullfile(cfg.leGUIRoot, 'tpm', 'labels_Neuromorphometrics.nii');
if ~isfile(petFile)
    fprintf('[PHG] Receptor surface figure skipped: %s missing.\n', petFile);
    return
end
if ~isfile(atlasFile)
    fprintf('[PHG] Receptor surface figure skipped: %s missing.\n', atlasFile);
    return
end

surface = phg.buildTransparentBrainSurface(cfg);

%% Sample the PET map at every cortical vertex
petInfo = niftiinfo(petFile);
petVolume = double(niftiread(petFile));
vertexDensity = localSample(petVolume, petInfo, surface.vertices);
% Robust limits: a few voxels near the sinuses carry extreme values that would
% otherwise flatten the whole cortical range.
limits = prctile(vertexDensity(isfinite(vertexDensity)), [2 98]);

%% Subcortical meshes from the atlas LeGUI already ships
atlasInfo = niftiinfo(atlasFile);
atlasVolume = niftiread(atlasFile);
structures = struct( ...
    'name',  {'Hippocampus', 'Amygdala'}, ...
    'index', {[47 48], [31 32]}, ...
    'color', {[0.06 0.42 0.47], [0.72 0.45 0.20]});

meshes = cell(numel(structures), 1);
for k = 1:numel(structures)
    mask = ismember(double(atlasVolume), structures(k).index);
    mask = smooth3(double(mask), 'gaussian', 5, 1.1);
    [faces, vertices] = isosurface(mask, 0.45);
    if isempty(faces)
        continue
    end
    % isosurface returns columns as (x,y,z) = (dim2,dim1,dim3); swap back
    % before applying the voxel-to-world transform.
    vertices = vertices(:, [2 1 3]);
    affine = atlasInfo.Transform.T';
    world = affine * [(vertices - 1)'; ones(1, size(vertices, 1))];
    meshes{k} = struct('faces', faces, 'vertices', world(1:3, :)', ...
        'color', structures(k).color, 'name', structures(k).name);
end
meshes = meshes(~cellfun(@isempty, meshes));

%% Contacts, grouped into leads
coordinates = channel.XYZMNI;
hasExcursion = channel.RespAreaNet ~= 0;
isDilation = channel.RespAreaNet > 0;
shaft = string(channel.PtID) + "|" + phg.parseLeadLabel(string(channel.Label));
isHippocampal = phg.cleanRegionLabels(channel.NMM) == "Hippocampus";

fig = figure('Color', 'w', 'Units', 'centimeters', 'Visible', 'off', ...
    'Position', [2 2 18.0 7.4]);
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
phg.setSafeLayout(layout);
% The panel labels sit above the axes box, so the layout gives up a little
% height at the top; without it they touch the export boundary.
layout.OuterPosition = [0.015 0.005 0.980 0.930];

viewSpec = struct('camera', [-1 0 0], 'up', [0 0 1]);

for panel = 1:2
    ax = nexttile(layout, panel);
    hold(ax, 'on');

    % Panel a is the map, so its surface is opaque. Panel b is the anatomy, so
    % its surface is a glass shell: at 0.92 alpha the subcortical meshes and
    % the depth contacts are simply invisible, which defeats the figure.
    if panel == 1
        patch(ax, 'Faces', surface.faces, 'Vertices', surface.vertices, ...
            'FaceVertexCData', vertexDensity, 'FaceColor', 'interp', ...
            'FaceAlpha', 0.97, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
            'AmbientStrength', 0.35, 'DiffuseStrength', 0.75, ...
            'SpecularStrength', 0.08, 'SpecularExponent', 12);
        colormap(ax, phg.divergingMap(256));
        clim(ax, limits);
    else
        patch(ax, 'Faces', surface.faces, 'Vertices', surface.vertices, ...
            'FaceColor', [0.62 0.64 0.67], 'FaceAlpha', 0.10, ...
            'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
            'AmbientStrength', 0.55, 'DiffuseStrength', 0.45, ...
            'SpecularStrength', 0.04);
    end

    for k = 1:numel(meshes) * (panel == 2)
        patch(ax, 'Faces', meshes{k}.faces, 'Vertices', meshes{k}.vertices, ...
            'FaceColor', meshes{k}.color, 'FaceAlpha', 1.0, ...
            'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
            'AmbientStrength', 0.45, 'DiffuseStrength', 0.85, ...
            'SpecularStrength', 0.25, 'SpecularExponent', 20);
    end

    if panel == 2
        localDrawLeads(ax, coordinates, shaft, hasExcursion, isDilation, isHippocampal, style);
    end

    daspect(ax, [1 1 1]);
    axis(ax, 'vis3d', 'off');
    view(ax, viewSpec.camera);
    camup(ax, viewSpec.up);
    camproj(ax, 'orthographic');
    camlight(ax, 'headlight');
    lighting(ax, 'gouraud');
    material(ax, 'dull');

    % Axis limits are pinned to the SURFACE, not to everything drawn. Some
    % contacts sit outside the template, and under `axis vis3d` they would
    % expand the bounding sphere the camera fits — which is why panel b
    % rendered smaller than panel a when it alone carried the contacts.
    bounds = [min(surface.vertices, [], 1); max(surface.vertices, [], 1)];
    xlim(ax, bounds(:,1)'); ylim(ax, bounds(:,2)'); zlim(ax, bounds(:,3)');
    span = diff(bounds, 1, 1);
    camzoom(ax, 0.62 * norm(span) / max(span([2 3])));
    % Both panels must render at one scale. camzoom is applied per axes and
    % the two axes do not end up with identical camera view angles, so the
    % first panel's angle is copied to the second rather than trusted to match.
    if panel == 1
        sharedViewAngle = ax.CameraViewAngle;
    else
        ax.CameraViewAngle = sharedViewAngle;
    end

    if panel == 1
        title(ax, 'Norepinephrine transporter density');
        ax.Tag = 'mapPanel';
        phg.addPanelLabel(ax, 'a');
    else
        title(ax, 'Contacts by coupling direction');
        phg.addPanelLabel(ax, 'b');
    end
    ax.Title.FontName = style.fontName;
    ax.Title.FontSize = style.titleFontSize;
    ax.Title.FontWeight = 'normal';
end

bar = colorbar(findobj(layout, 'Type', 'axes', 'Tag', 'mapPanel'));
bar.Layout.Tile = 'south';
bar.Label.String = 'NET density (atlas units)';
bar.Label.FontName = style.fontName;
bar.Label.FontSize = style.labelFontSize;
bar.FontName = style.fontName;
bar.FontSize = style.tickFontSize;
bar.Ticks = limits;
bar.TickLabels = compose('%.2f', limits);

qc = phg.exportPublicationFigure(fig, cfg, "Fig10_receptor_surface", ...
    'pdfContentType', "image");
close(fig);
end

% -------------------------------------------------------------------------
function values = localSample(volume, info, coordinates)
affine = info.Transform.T';
voxel = round(affine \ [coordinates'; ones(1, size(coordinates, 1))])';
voxel = voxel(:, 1:3) + 1;
inside = all(voxel >= 1, 2) & voxel(:,1) <= info.ImageSize(1) & ...
    voxel(:,2) <= info.ImageSize(2) & voxel(:,3) <= info.ImageSize(3);
values = nan(size(coordinates, 1), 1);
values(inside) = volume(sub2ind(info.ImageSize, voxel(inside,1), ...
    voxel(inside,2), voxel(inside,3)));
end

function localDrawLeads(ax, coordinates, shaft, hasExcursion, isDilation, isHippocampal, style)
%LOCALDRAWLEADS Join contacts on one shaft, then mark each contact.
% Drawing all 146 shafts overlaid in one template space is visual noise, not
% information. Only shafts that reach the hippocampus are traced, which is
% the sampling the paper's primary contrast depends on.
mesialShafts = unique(shaft(isHippocampal));
for k = 1:numel(mesialShafts)
    on = shaft == mesialShafts(k);
    if sum(on) < 2
        continue
    end
    points = coordinates(on, :);
    % Order along the shaft's principal axis so the connecting line follows
    % the trajectory rather than the arbitrary label order.
    centred = points - mean(points, 1);
    [~, ~, basis] = svd(centred, 'econ');
    [~, order] = sort(centred * basis(:,1));
    points = points(order, :);
    plot3(ax, points(:,1), points(:,2), points(:,3), '-', ...
        'Color', [0.35 0.35 0.38 0.55], 'LineWidth', 1.1);
end

quiet = ~hasExcursion;
scatter3(ax, coordinates(quiet,1), coordinates(quiet,2), coordinates(quiet,3), ...
    8, style.lightGray, 'filled', 'MarkerFaceAlpha', 0.22);
dilating = hasExcursion & isDilation;
scatter3(ax, coordinates(dilating,1), coordinates(dilating,2), ...
    coordinates(dilating,3), 26, style.dilation, 'filled', ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.3);
constricting = hasExcursion & ~isDilation;
scatter3(ax, coordinates(constricting,1), coordinates(constricting,2), ...
    coordinates(constricting,3), 26, style.constriction, 'filled', ...
    'MarkerEdgeColor', 'white', 'LineWidth', 0.3);
end
