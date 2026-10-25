function [faces, vertices] = mesh_to_triangles(mesh, fib)
%MESH_TO_TRIANGLES Convert an edge mesh into a triangulated surface of cylinders.
%   [faces, vertices] = mesh_to_triangles(mesh, fib)
%
%   Every edge of the mesh is replaced by a cylinder of radius fib.cylRadius
%   (optionally closed by hemispherical caps, see cylinder_between). The
%   cylinders are built independently and concatenated: vertices are not
%   shared between different edges, and intersecting cylinders are not merged.
%   Face normals point outwards from each cylinder.
%
%   Inputs
%     mesh - struct with fields
%              vertices : Mx3 coordinates of the edge end points
%              edges    : Nx2 vertex indices [start end] of each edge
%     fib  - fibre discretisation struct (from params_reference) with fields
%              cylRadius   : cylinder radius
%              nang        : number of angular segments around the axis
%              ringSpacing : maximum axial distance between rings; edge i
%                            gets max(2, ceil(L_i/ringSpacing) + 1) rings
%              capEnds     : true to close each cylinder with hemispheres
%              nCapRings   : number of latitude rings per hemisphere
%                            (used only if capEnds is true)
%
%   Outputs
%     faces    - Fx3 vertex indices (one triangle per row)
%     vertices - Vx3 vertex coordinates
%
%   Typical usage (see build_fibres_aligned.m):
%     mesh  = mesh_fibres_aligned(cfg, angle);
%     [F,V] = mesh_to_triangles(mesh, fib);
%     stlwrite(triangulation(F, V), 'parallel_fibres.stl');

% No edges: return an empty triangulation.
if isempty(mesh.edges)
    faces = zeros(0,3);
    vertices = zeros(0,3);
    return
end

edges = mesh.edges;
verts = mesh.vertices;

if size(edges,2) ~= 2 || size(verts,2) ~= 3
    error('Invalid mesh format: edges must be Nx2 and vertices must be Mx3.');
end

% Exact total size, to avoid reallocations (same counts as in
% cylinder_between, including its rounding of nang). Edges of zero length
% produce no faces; the arrays are trimmed after the loop.
% Number of rings per edge: enough for a spacing <= fib.ringSpacing.
edgeLengths = sqrt(sum((verts(edges(:,2),:) - verts(edges(:,1),:)).^2, 2));
nzPerEdge = max(2, ceil(edgeLengths ./ fib.ringSpacing) + 1);

nang = max(3, round(fib.nang));
if fib.capEnds
    capFaces = 2*nang*(2*fib.nCapRings - 1);   % two hemispheres
    capVerts = 2*(1 + (fib.nCapRings-1)*nang); % pole + inner rings (rim shared), two hemispheres
else
    capFaces = 0;
    capVerts = 0;
end

facesPerEdge = 2*nang*(nzPerEdge-1) + capFaces;
vertsPerEdge = nang*nzPerEdge + capVerts;

faces = zeros(sum(facesPerEdge), 3);
vertices = zeros(sum(vertsPerEdge), 3);

% Build each cylinder in turn and append it. Its local vertex indices are
% shifted by the number of vertices already stored (vCount), since each
% cylinder has its own vertices.
fCount = 0;   % faces stored so far
vCount = 0;   % vertices stored so far

for i = 1:size(edges,1)
    p1 = verts(edges(i,1),:);
    p2 = verts(edges(i,2),:);

    nzEdge = nzPerEdge(i);

    [fLocal, vLocal] = cylinder_between(p1, p2, fib.cylRadius, fib.nang, nzEdge, fib.capEnds, fib.nCapRings);

    % Zero-length edge (or r <= 0): nothing to add.
    if isempty(fLocal)
        continue
    end

    nv = size(vLocal,1);
    nf = size(fLocal,1);

    vertices(vCount + (1:nv), :) = vLocal;
    faces(fCount + (1:nf), :) = fLocal + vCount;

    vCount = vCount + nv;
    fCount = fCount + nf;
end

% Remove the rows reserved for zero-length edges, if any.
vertices = vertices(1:vCount,:);
faces = faces(1:fCount,:);
end
