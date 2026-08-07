function [faces, vertices] = leadCylinder(startPoint, endPoint, radius, options)
%LEADCYLINDER Closed cylinder between two points, as a patch.
%
%   Reimplements the geometry of Lead-DBS's ea_cylinder for a patch surface
%   rather than a mesh object, so that a shaft or a contact composes with the
%   patch-based rendering used everywhere else in this project.
%
%   Lead-DBS builds its cylinder along one axis and rotates it onto the
%   required direction. The same is done here with an explicit orthonormal
%   frame instead of MATLAB's `rotate`, which acts on a graphics handle and
%   cannot be used before the object exists.
%
%   [F, V] = phg.leadCylinder(P1, P2, R) returns faces and vertices of a
%   closed cylinder of radius R whose axis runs from P1 to P2.

arguments
    startPoint (1,3) double
    endPoint (1,3) double
    radius (1,1) double {mustBePositive}
    options.Facets (1,1) double {mustBeInteger, mustBePositive} = 24
end

axisVector = endPoint - startPoint;
axisLength = norm(axisVector);
if axisLength == 0
    faces = zeros(0, 3);
    vertices = zeros(0, 3);
    return
end
unitAxis = axisVector ./ axisLength;

% Any vector not parallel to the axis gives a starting point for the frame;
% picking the world axis the shaft is *least* aligned with keeps the cross
% product well conditioned for shafts that run along a cardinal direction.
[~, weakest] = min(abs(unitAxis));
seed = zeros(1, 3);
seed(weakest) = 1;
firstNormal = cross(unitAxis, seed);
firstNormal = firstNormal ./ norm(firstNormal);
secondNormal = cross(unitAxis, firstNormal);

angle = linspace(0, 2*pi, options.Facets + 1)';
angle(end) = [];
ring = radius .* (cos(angle) .* firstNormal + sin(angle) .* secondNormal);

bottomRing = startPoint + ring;
topRing = endPoint + ring;
vertices = [bottomRing; topRing; startPoint; endPoint];

n = options.Facets;
next = [2:n, 1];
% Side wall as two triangles per facet, then a triangle fan for each cap.
sideFaces = [ (1:n)', next', (next + n)'; (1:n)', (next + n)', ((1:n) + n)' ];
bottomCap = [ repmat(2*n + 1, n, 1), next', (1:n)' ];
topCap = [ repmat(2*n + 2, n, 1), ((1:n) + n)', (next + n)' ];
faces = [sideFaces; bottomCap; topCap];
end
