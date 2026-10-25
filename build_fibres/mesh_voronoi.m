function [mesh, info] = mesh_voronoi(points, limits, enableBoundaryCollinearCollapse, mergeVertexTolerance)
%VORONOICUBEMESH_V2 Bounded 3-D Voronoi mesh of seeds in a box.
%   mesh = voronoiCubeMesh_v2(points)
%   mesh = voronoiCubeMesh_v2(points, limits)
%   mesh = voronoiCubeMesh_v2(points, limits, enableBoundaryCollinearCollapse)
%   mesh = voronoiCubeMesh_v2(points, limits, enableBoundaryCollinearCollapse, mergeVertexTolerance)
%   [mesh, info] = voronoiCubeMesh_v2(...)
%
%   Same inputs and output fields as voronoiCubeMesh, with the issues of
%   that version fixed (see notes/triangulation_geometry.md, section 5).
%   mesh.vertices and mesh.edges can be passed to mesh_to_triangles.
%
%   INPUT
%   points : Nx3 array of distinct seed coordinates, strictly inside the box.
%   limits : Optional struct defining the box with fields
%            xmin, xmax, ymin, ymax, zmin, zmax.
%            If omitted or empty, defaults to the unit cube [0,1]^3.
%   enableBoundaryCollinearCollapse : Optional logical flag (default false).
%            If true, remove degree-2 vertices of mesh.edges that lie on the
%            box surface and whose two edges are within 5 degrees of
%            collinear, joining their two neighbours. With the
%            construction used here such vertices are rare, so this usually
%            changes nothing.
%   mergeVertexTolerance : Optional nonnegative scalar (default 0 =
%            disabled), an ABSOLUTE length in the units of limits (e.g. 0.1
%            is 10% of the unit cube but 0.2% of a 50-unit box). If > 0,
%            vertices of mesh.edges connected by chains of pairwise
%            distances <= mergeVertexTolerance are merged into one vertex at
%            the centroid of the chain (kept on the box faces it touches).
%            Useful to remove edges much shorter than the cylinder diameter
%            used later in mesh_to_triangles.
%
%   METHOD
%   1) Reflect the seeds across the six box faces and compute the Voronoi
%      diagram of the 7N points (voronoin). For x in the box and a seed p,
%      the reflection p' of p across any face satisfies |x - p'| >= |x - p|,
%      so the cell of an original seed is exactly its Voronoi cell clipped
%      to the box: no further clipping is needed.
%   2) Merge the near-coincident vertices produced by Qhull (relative
%      tolerance), so that each geometric vertex has one id.
%   3) Two cells sharing >= 3 vertices share a face (the shared vertices
%      lie on the bisector plane of their seeds). A face between seeds i
%      and j <= N is internal; a face between seed i and a reflected seed
%      lies on the box surface. Each face loop is ordered by angle around
%      its centroid, anticlockwise about the normal (p_j - p_i)/|p_j - p_i|
%      (i.e. outward from cell i).
%   4) Edges are the consecutive vertex pairs of the face loops.
%   5) Optional post-processing (collinear collapse, vertex merging); the
%      faces, cells and edge lists are updated consistently.
%
%   OUTPUT (mesh struct)
%   mesh.points      : Input seeds (Nx3).
%   mesh.vertices    : Vertices used by mesh.edges (Mused x 3).
%   mesh.verticesAll : All vertices of the cells (Mall x 3), after the
%                      optional post-processing; mesh.cells, mesh.cellFaces,
%                      mesh.faces, mesh.facesAll and mesh.edgesAll index
%                      into this list.
%   mesh.cells       : Cell array, one entry per seed: vertex ids of its cell.
%   mesh.cellFaces   : Cell array, one entry per seed: its internal faces
%                      (cell array of vertex-id loops).
%   mesh.faces       : Internal faces (shared by two cells), as a cell array
%                      of row vectors of vertex ids.
%   mesh.edges       : Edges of the internal faces, indexing mesh.vertices.
%                      They include the rims where the internal faces meet
%                      the box surface (these lie ON the box faces, so
%                      cylinders built on them by mesh_to_triangles extend
%                      outside the box by their radius; to avoid this, pass
%                      limits reduced by the cylinder radius). Segments
%                      lying on the 12 box edges are excluded.
%   mesh.facesAll    : Internal faces followed by the faces on the box surface.
%   mesh.edgesAll    : Edges of mesh.facesAll (including box-edge segments),
%                      indexing mesh.verticesAll.
%
%   info (optional second output) : struct with diagnostics
%   info.nMissingCells   : seeds without a valid cell (should be 0)
%   info.closureError    : max over cells of |sum of outward face area
%                          vectors| / L^2 (0 for closed cells; L = largest
%                          box side)
%   info.volumeError     : |volume enclosed by the box-surface faces - box
%                          volume| / box volume (0 if the box is covered)
%   info.minEdgeLength   : shortest edge of mesh.edges
%   info.nCollapsed      : vertices removed by the collinear collapse
%   info.nMerged         : vertices removed by the merging

%% Inputs and defaults
if nargin < 1 || isempty(points) || size(points,2) ~= 3
    error('voronoiCubeMesh_v2:InvalidPoints', 'points must be a non-empty Nx3 array.');
end
if ~all(isfinite(points(:)))
    error('voronoiCubeMesh_v2:InvalidPoints', 'points must be finite.');
end

if nargin < 2 || isempty(limits)
    limits.xmin = 0;
    limits.xmax = 1;
    limits.ymin = 0;
    limits.ymax = 1;
    limits.zmin = 0;
    limits.zmax = 1;
end

if nargin < 3 || isempty(enableBoundaryCollinearCollapse)
    enableBoundaryCollinearCollapse = false;
end

if nargin < 4 || isempty(mergeVertexTolerance)
    mergeVertexTolerance = 0;
end

if ~isscalar(mergeVertexTolerance) || mergeVertexTolerance < 0
    error('voronoiCubeMesh_v2:InvalidMergeTolerance', ...
        'mergeVertexTolerance must be a nonnegative scalar.');
end

requiredFields = {'xmin','xmax','ymin','ymax','zmin','zmax'};
for f = 1:numel(requiredFields)
    if ~isfield(limits,requiredFields{f})
        error('voronoiCubeMesh_v2:InvalidLimits','Missing limits.%s field.',requiredFields{f});
    end
end

if limits.xmin >= limits.xmax || limits.ymin >= limits.ymax || limits.zmin >= limits.zmax
    error('voronoiCubeMesh_v2:InvalidLimits','Each limits min value must be strictly smaller than its max value.');
end

n  = size(points,1);
lo = [limits.xmin limits.ymin limits.zmin];
hi = [limits.xmax limits.ymax limits.zmax];

% Tolerances relative to the box size L (the old version used absolute ones).
L         = max(hi - lo);
tolVertex = 1e-9  * L;      % merging of near-coincident Qhull vertices / seeds
tolPlane  = 1e-9  * L;      % "point lies on a box plane"
tolArea   = 1e-12 * L^2;    % faces with smaller area are discarded

% The reflection argument requires every seed strictly inside the box: a
% seed outside leaves part of the box to a reflected seed (a hole), and a
% seed on a face coincides with its own reflection.
distToBoundary = min([points - lo, hi - points], [], 2);
nBad = nnz(distToBoundary <= tolPlane);
if nBad > 0
    error('voronoiCubeMesh_v2:SeedsNotInside', ...
        '%d seed(s) lie outside or on the boundary of the box; all seeds must be strictly inside.', nBad);
end

% Coincident seeds would give a degenerate Voronoi diagram.
if size(uniquetol(points, tolVertex, 'ByRows', true, 'DataScale', 1), 1) < n
    error('voronoiCubeMesh_v2:DuplicateSeeds', 'points contains coincident seeds.');
end

%% 1) Voronoi diagram of the seeds and their reflections
% Block k = 1..6 holds the reflections across xmin, xmax, ymin, ymax,
% zmin, zmax: the reflection of seed i across face k is row k*n + i.
p_ext = [
    points;
    [2*limits.xmin-points(:,1) points(:,2) points(:,3)];
    [2*limits.xmax-points(:,1) points(:,2) points(:,3)];
    [points(:,1) 2*limits.ymin-points(:,2) points(:,3)];
    [points(:,1) 2*limits.ymax-points(:,2) points(:,3)];
    [points(:,1) points(:,2) 2*limits.zmin-points(:,3)];
    [points(:,1) points(:,2) 2*limits.zmax-points(:,3)];
];
nExt = size(p_ext,1);

[V_ext, C_ext] = voronoin(p_ext);

%% 2) Merge near-coincident vertices
% V_ext(1,:) is the point at infinity (used only by unbounded cells of
% reflected seeds); it gets id 0 and is dropped.
[U, ~, ic] = uniquetol(V_ext(2:end,:), tolVertex, 'ByRows', true, 'DataScale', 1);
vid = [0; ic(:)];   % vid(k) = id in U of Voronoi vertex k

cellV = cell(nExt,1);   % merged vertex ids of each cell (column, sorted)
for c = 1:nExt
    g = vid(C_ext{c});
    cellV{c} = unique(g(g > 0));
    cellV{c} = cellV{c}(:);
end

%% 3) Faces from pairs of cells sharing >= 3 vertices
% Sparse incidence matrix (cell x vertex); S(a,b) = number of vertices
% shared by original cell a and cell b.
nPerCell = cellfun(@numel, cellV);
M = sparse(repelem((1:nExt).', nPerCell), vertcat(cellV{:}), 1, nExt, size(U,1));
S = M(1:n,:) * M.';
[aList, bList, nShared] = find(S);

% Keep each internal pair once (b > a); keep all pairs with a reflection.
keep  = nShared >= 3 & (bList > n | bList > aList);
aList = aList(keep);
bList = bList(keep);

nCand = numel(aList);
loops = cell(nCand,1);
areaVec = zeros(nCand,3);   % area vector of each face, oriented outward from cell a
valid = false(nCand,1);

for f = 1:nCand
    a = aList(f);
    b = bList(f);
    ids = intersect(cellV{a}, cellV{b});
    nrm = p_ext(b,:) - p_ext(a,:);
    nrm = nrm / norm(nrm);

    [loop, area] = orderFaceLoop(ids, U, nrm);
    if area > tolArea
        valid(f) = true;
        loops{f} = loop;
        areaVec(f,:) = area * nrm;
    end
end

loops   = loops(valid);
aList   = aList(valid);
bList   = bList(valid);
areaVec = areaVec(valid,:);
isBoundaryFace = bList > n;

%% Diagnostics: closed cells and full coverage of the box
% A closed surface has zero total area vector, so each cell's outward face
% area vectors (+ for cell a, - for cell b) must sum to zero. The volume
% enclosed by the box-surface faces, (1/3) * sum(x . n A) over them, must be
% the box volume.
nMissingCells = nnz(cellfun(@numel, cellV(1:n)) < 4);
sumArea = zeros(n,3);
for d = 1:3
    sumArea(:,d) = accumarray(aList, areaVec(:,d), [n 1]);
    inner = ~isBoundaryFace;
    sumArea(:,d) = sumArea(:,d) - accumarray(bList(inner), areaVec(inner,d), [n 1]);
end
closureError = max(vecnorm(sumArea, 2, 2)) / L^2;

boundaryLoops = loops(isBoundaryFace);
boundaryArea  = areaVec(isBoundaryFace,:);
volBoundary = 0;
for f = 1:numel(boundaryLoops)
    volBoundary = volBoundary + dot(U(boundaryLoops{f}(1),:), boundaryArea(f,:)) / 3;
end
volBox = prod(hi - lo);
volumeError = abs(volBoundary - volBox) / volBox;

if nMissingCells > 0
    warning('voronoiCubeMesh_v2:MissingCells', '%d seed(s) have no valid cell.', nMissingCells);
end
if closureError > 1e-8 || volumeError > 1e-8
    warning('voronoiCubeMesh_v2:InconsistentCells', ...
        'Cells not closed or box not covered (closure error %.2e, volume error %.2e).', ...
        closureError, volumeError);
end

%% Global vertex list: vertices of the original cells only
usedU = unique(vertcat(cellV{1:n}));
uToAll = zeros(size(U,1),1);
uToAll(usedU) = 1:numel(usedU);
verticesAll = U(usedU,:);

loops = cellfun(@(lp) reshape(uToAll(lp), 1, []), loops, 'UniformOutput', false);
cells = cellfun(@(cv) uToAll(cv), cellV(1:n), 'UniformOutput', false);

%% 4) Edges of the internal faces, without box-edge segments
edges = removeCubeEdgeSegments(edgesFromLoops(loops(~isBoundaryFace)), verticesAll, lo, hi, tolPlane);

%% 5) Optional post-processing (faces and cells updated with the same maps)
nCollapsed = 0;
if enableBoundaryCollinearCollapse
    [edges, removed] = collapseNearCollinearDegree2Vertices(edges, verticesAll, lo, hi, tolPlane, 5);
    nCollapsed = nnz(removed);
    if nCollapsed > 0
        loops = cellfun(@(lp) lp(~removed(lp)), loops, 'UniformOutput', false);
        cells = cellfun(@(cv) cv(~removed(cv)), cells, 'UniformOutput', false);
    end
end

nMerged = 0;
if mergeVertexTolerance > 0
    [edges, verticesAll, oldToMerged] = mergeCloseGraphVertices(edges, verticesAll, ...
        mergeVertexTolerance, lo, hi, tolPlane);
    nMerged = nnz(oldToMerged ~= (1:numel(oldToMerged)).');
    if nMerged > 0
        loops = cellfun(@(lp) oldToMerged(lp(:)).', loops, 'UniformOutput', false);
        cells = cellfun(@(cv) unique(oldToMerged(cv)), cells, 'UniformOutput', false);
    end
end

if nCollapsed > 0 || nMerged > 0
    % Merged vertices snapped onto two box planes can create new segments
    % along a box edge: remove them as in step 4.
    edges = removeCubeEdgeSegments(edges, verticesAll, lo, hi, tolPlane);

    % Remove repeated consecutive ids (cyclically) and faces left with
    % fewer than 3 distinct vertices.
    for f = 1:numel(loops)
        lp = loops{f};
        if ~isempty(lp)
            lp = lp(lp ~= lp([end 1:end-1]));
        end
        if numel(unique(lp)) < 3
            lp = [];
        end
        loops{f} = lp;
    end
    ok = ~cellfun(@isempty, loops);
    loops = loops(ok);
    aList = aList(ok);
    bList = bList(ok);
    isBoundaryFace = isBoundaryFace(ok);
end

%% Assemble outputs
faces    = loops(~isBoundaryFace);
facesAll = [faces; loops(isBoundaryFace)];
edgesAll = edgesFromLoops(facesAll);

% Internal faces of each cell (each internal face belongs to cells a and b).
% Sort the (cell, face) pairs by cell and split them into one block per cell.
aInt = aList(~isBoundaryFace);
bInt = bList(~isBoundaryFace);
faceIdx = (1:numel(faces)).';
owner = sortrows([[aInt; bInt] [faceIdx; faceIdx]]);
fidx = owner(:,2);
counts = accumarray(owner(:,1), 1, [n 1]);
cellFaceIdx = mat2cell(fidx, counts, 1);
cellFaces = cellfun(@(ix) faces(ix), cellFaceIdx, 'UniformOutput', false);

% Compact list of the vertices used by mesh.edges.
if isempty(edges)
    vertices = zeros(0,3);
    edges = zeros(0,2);
else
    usedVertexIds = unique(edges(:));
    oldToNew = zeros(size(verticesAll,1),1);
    oldToNew(usedVertexIds) = 1:numel(usedVertexIds);
    vertices = verticesAll(usedVertexIds,:);
    edges = reshape(oldToNew(edges), [], 2);   % reshape: keeps Nx2 also when N = 1
end

mesh.points = points;
mesh.vertices = vertices;
mesh.verticesAll = verticesAll;
mesh.cells = cells;
mesh.cellFaces = cellFaces;
mesh.faces = faces;
mesh.edges = edges;
mesh.facesAll = facesAll;
mesh.edgesAll = edgesAll;

if isempty(edges)
    minEdgeLength = NaN;
else
    minEdgeLength = min(vecnorm(vertices(edges(:,2),:) - vertices(edges(:,1),:), 2, 2));
end
info = struct('nMissingCells', nMissingCells, 'closureError', closureError, ...
    'volumeError', volumeError, 'minEdgeLength', minEdgeLength, ...
    'nCollapsed', nCollapsed, 'nMerged', nMerged);

end

function [loop, area] = orderFaceLoop(ids, vertices, nrm)
% Order the vertices ids of a planar convex face with unit normal nrm by
% angle around their centroid, anticlockwise about nrm. Also return the
% polygon area (shoelace formula in the plane of the face).

pts = vertices(ids,:);
c = mean(pts,1);

% Orthonormal basis (u, v) of the face plane with cross(u, v) = nrm.
if abs(nrm(1)) < 0.9
    ref = [1 0 0];
else
    ref = [0 1 0];
end
u = cross(nrm, ref);
u = u / norm(u);
v = cross(nrm, u);

rel = pts - c;
q = [rel*u.' rel*v.'];
[~, order] = sort(atan2(q(:,2), q(:,1)));

loop = reshape(ids(order), 1, []);
q = q(order,:);
area = 0.5 * sum(q(:,1).*q([2:end 1],2) - q([2:end 1],1).*q(:,2));

end

function edges = edgesFromLoops(loops)
% Unique undirected edges (sorted rows) from closed vertex loops.

if isempty(loops)
    edges = zeros(0,2);
    return
end

pairs = cellfun(@(lp) [lp(:) lp([2:end 1]).'], loops, 'UniformOutput', false);
edges = unique(sort(vertcat(pairs{:}), 2), 'rows');

end

function edges = removeCubeEdgeSegments(edges, vertices, lo, hi, tol)
% Remove edges whose two end points both lie on the same box edge (i.e. on
% the same pair of box planes).

if isempty(edges)
    return
end

p1 = vertices(edges(:,1),:);
p2 = vertices(edges(:,2),:);

% onPlane(:, d, 1) / (:, d, 2): both end points on the min / max plane of dimension d.
onPlane = false(size(edges,1), 3, 2);
for d = 1:3
    onPlane(:,d,1) = abs(p1(:,d) - lo(d)) <= tol & abs(p2(:,d) - lo(d)) <= tol;
    onPlane(:,d,2) = abs(p1(:,d) - hi(d)) <= tol & abs(p2(:,d) - hi(d)) <= tol;
end

% On a box edge = on a plane in at least two different dimensions.
onCubeEdge = sum(any(onPlane, 3), 2) >= 2;
edges = edges(~onCubeEdge,:);

end

function tf = isOnBoundarySurface(p, lo, hi, tol)
% True for the rows of p that lie on at least one box plane.

tf = any(abs(p - lo) <= tol | abs(p - hi) <= tol, 2);

end

function [edges, removed] = collapseNearCollinearDegree2Vertices(edges, vertices, lo, hi, tol, angleTolDeg)
% Remove degree-2 vertices on the box surface whose two edges make an angle
% within angleTolDeg of 180 degrees, joining their two neighbours. Uses
% adjacency lists and a work queue: when a vertex is removed, its two
% neighbours are re-examined, since their angles have changed.

nV = size(vertices,1);
removed = false(nV,1);
if isempty(edges)
    return
end

nbr = cell(nV,1);
for e = 1:size(edges,1)
    nbr{edges(e,1)}(end+1) = edges(e,2);
    nbr{edges(e,2)}(end+1) = edges(e,1);
end

onB = isOnBoundarySurface(vertices, lo, hi, tol);
cosThreshold = -cosd(angleTolDeg);

queue = find(cellfun(@numel, nbr) == 2 & onB);
k = 1;
while k <= numel(queue)
    v = queue(k);
    k = k + 1;

    if removed(v) || numel(nbr{v}) ~= 2 || ~onB(v)
        continue
    end

    n1 = nbr{v}(1);
    n2 = nbr{v}(2);
    a = vertices(n1,:) - vertices(v,:);
    b = vertices(n2,:) - vertices(v,:);
    na = norm(a);
    nb = norm(b);
    if na <= tol || nb <= tol || dot(a,b)/(na*nb) > cosThreshold
        continue
    end

    % Replace n1 - v - n2 by n1 - n2 (unless that edge already exists).
    nbr{n1}(nbr{n1} == v) = [];
    nbr{n2}(nbr{n2} == v) = [];
    if ~any(nbr{n1} == n2)
        nbr{n1}(end+1) = n2;
        nbr{n2}(end+1) = n1;
    end
    nbr{v} = [];
    removed(v) = true;
    queue = [queue; n1; n2]; %#ok<AGROW>
end

% Rebuild the edge list from the adjacency lists (each edge once).
nDeg = cellfun(@numel, nbr);
I = repelem((1:nV).', nDeg);
J = [nbr{:}].';
edges = sortrows([I(I < J) J(I < J)]);

end

function [edges, vertices, oldToMerged] = mergeCloseGraphVertices(edges, vertices, mergeTol, lo, hi, tolPlane)
% Merge the vertices of edges connected by chains of pairwise distances
% <= mergeTol. Each group is moved to its centroid, snapped back onto any
% box plane touched by one of its members; all ids of a group are mapped to
% its smallest id (oldToMerged) and the other ids become unused.

nV = size(vertices,1);
oldToMerged = (1:nV).';
if isempty(edges)
    return
end

used = unique(edges(:));
pts = vertices(used,:);
nUsed = numel(used);

% Candidate pairs by a sweep in x: after sorting by x, only points with
% x-distance <= mergeTol need to be compared.
[xs, order] = sort(pts(:,1));
parent = (1:nUsed).';
for ii = 1:nUsed-1
    i = order(ii);
    jj = ii + 1;
    while jj <= nUsed && xs(jj) - xs(ii) <= mergeTol
        j = order(jj);
        if norm(pts(i,:) - pts(j,:)) <= mergeTol
            [parent, ri] = ufFind(parent, i);
            [parent, rj] = ufFind(parent, j);
            if ri ~= rj
                parent(max(ri,rj)) = min(ri,rj);
            end
        end
        jj = jj + 1;
    end
end

roots = zeros(nUsed,1);
for i = 1:nUsed
    [parent, roots(i)] = ufFind(parent, i);
end

% Groups of local indices sharing a root.
[~, ~, g] = unique(roots);
groups = accumarray(g, (1:nUsed).', [], @(x) {sort(x)});

for k = 1:numel(groups)
    members = used(groups{k});   % global ids, increasing
    if numel(members) < 2
        continue
    end

    newPos = mean(vertices(members,:), 1);
    P = vertices(members,:);
    for d = 1:3
        onLo = any(abs(P(:,d) - lo(d)) <= tolPlane);
        onHi = any(abs(P(:,d) - hi(d)) <= tolPlane);
        if onLo && ~onHi
            newPos(d) = lo(d);
        elseif onHi && ~onLo
            newPos(d) = hi(d);
        end
    end

    vertices(members,:) = repmat(newPos, numel(members), 1);
    oldToMerged(members) = members(1);
end

edges = reshape(oldToMerged(edges), [], 2);
edges = edges(edges(:,1) ~= edges(:,2),:);
edges = unique(sort(edges, 2), 'rows');

end

function [parent, root] = ufFind(parent, node)
% Union-find root of node, with path compression.

root = node;
while parent(root) ~= root
    root = parent(root);
end
while parent(node) ~= root
    next = parent(node);
    parent(node) = root;
    node = next;
end

end
