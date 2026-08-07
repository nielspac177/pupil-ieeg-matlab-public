function outline = glassBrainOutline(surface, projection, options)
%GLASSBRAINOUTLINE Silhouette of an MNI surface projected onto a viewing plane.
%   The outline is derived from the template surface that
%   phg.buildTransparentBrainSurface produces, rather than from a hand-traced
%   path, so the silhouette and the plotted contacts are guaranteed to live in
%   the same coordinate frame.
%
%   projection is a struct with fields:
%     u, v      1x3 unit row vectors giving the horizontal and vertical
%               screen axes in MNI coordinates
%     mask      optional logical selector applied to the surface vertices
%               before the silhouette is computed (for hemisphere views)
%
%   outline has fields u, v (closed polygon) and midline (an inner contour
%   drawn from paramedian vertices to hint at internal structure).

arguments
    surface (1,1) struct
    projection (1,1) struct
    options.MaxVertices (1,1) double = 9000
    options.ShrinkFactor (1,1) double = 0.55
    options.MidlineHalfWidth (1,1) double = 9
end

outline = struct('u', [], 'v', [], 'midlineU', [], 'midlineV', []);
vertices = surface.vertices;
if isempty(vertices)
    return
end

if isfield(projection, 'mask') && ~isempty(projection.mask)
    vertices = vertices(projection.mask(vertices), :);
end
if isempty(vertices)
    return
end

if size(vertices, 1) > options.MaxVertices
    step = ceil(size(vertices, 1) / options.MaxVertices);
    vertices = vertices(1:step:end, :);
end

u = vertices * projection.u(:);
v = vertices * projection.v(:);
index = boundary(u, v, options.ShrinkFactor);
outline.u = u(index);
outline.v = v(index);

% Inner contour: the silhouette of paramedian tissue. On lateral and superior
% views this traces the corpus callosum / brainstem region and keeps the
% otherwise-empty interior from reading as a flat blob.
normal = cross(projection.u(:)', projection.v(:)');
distance = vertices * normal(:);
near = abs(distance - median(distance)) <= options.MidlineHalfWidth;
if sum(near) > 50
    innerIndex = boundary(u(near), v(near), 0.75);
    innerU = u(near);
    innerV = v(near);
    outline.midlineU = innerU(innerIndex);
    outline.midlineV = innerV(innerIndex);
end
end
