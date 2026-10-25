function [faces, vertices] = cylinder_between(p1, p2, r, nang, nz, capEnds, nCapRings)
%CYLINDER_BETWEEN Triangulated cylinder connecting p1 and p2.
%   [F,V] = cylinder_between(p1,p2,r,nang,nz,false)
%   [F,V] = cylinder_between(p1,p2,r,nang,nz,true,nCapRings)
%
%   p1, p2    : 1x3 endpoints of the axis.
%   r         : radius (empty output if r <= 0 or p1 == p2).
%   nang      : number of angular segments (rounded, >= 3).
%   nz        : number of axial rings, equally spaced from p1 to p2
%               (rounded, >= 2).
%   capEnds   : logical, add hemispherical end caps of radius r centred at
%               p1 and p2 (required, no default).
%   nCapRings : number of latitude rings of each hemisphere, excluding the
%               pole and including the rim. Required when capEnds is true
%               (unused, and may be omitted, when capEnds is false).
%
%   The surface is built in a local frame with the axis along +z from 0
%   to L = |p2 - p1|, then rotated onto p2 - p1 (Rodrigues' formula) and
%   translated to p1.
%   All faces are ordered anticlockwise when seen from outside, so the
%   normals given by the right-hand rule point outwards.
%   The rim ring of each hemisphere is shared with the first/last ring of
%   the side surface (no duplicated vertices), so the capped cylinder is a
%   closed, watertight surface.
%
%   Sizes: nang*nz + 2*(1 + (nCapRings-1)*nang) vertices and
%   2*nang*(nz-1) + 2*nang*(2*nCapRings - 1) faces (cap terms only if
%   capEnds is true).


nang = max(3, round(nang));
nz = max(2, round(nz));

axisVec = p2 - p1;
L = norm(axisVec);

if L < eps || r <= 0
    faces = zeros(0,3);
    vertices = zeros(0,3);
    return
end

% Angular positions of the nang vertices of each ring and axial positions
% of the nz rings.
theta = (0:nang-1)' * (2*pi/nang);
ct = cos(theta);
st = sin(theta);
zvals = linspace(0,L,nz).';

% Side vertices in local frame (axis aligned with +z); ring iz occupies
% rows (iz-1)*nang + (1:nang).
vertices = [
    r * repmat(ct,nz,1), ...
    r * repmat(st,nz,1), ...
    kron(zvals,ones(nang,1))
];

% Each hemisphere has nang*(2*nCapRings - 1) faces (see hemisphereMesh).
nSideFaces = 2 * nang * (nz-1);
if capEnds
    nCapFaces = 2 * nang * (2*nCapRings - 1);
else
    nCapFaces = 0;
end
faces = zeros(nSideFaces + nCapFaces,3);

% Side faces: with i0 -> i0n anticlockwise about +z and i0 -> i1 along +z,
% the order (i0, i1n, i1) gives the outward normal
% (i1n - i0) x (i1 - i0) ~ e_theta x e_z = e_r.

% Between rings iz and iz+1, each angular segment (vertex j -> j+1, with
% wrap-around nang -> 1 via the "n" index vectors) is split into two
% triangles.
frow = 1;   % next free row of faces
for iz = 1:nz-1
    base0 = (iz-1) * nang;   % offset of ring iz
    base1 = iz * nang;       % offset of ring iz+1

    i0 = base0 + (1:nang);
    i1 = base1 + (1:nang);
    i0n = base0 + [2:nang 1];
    i1n = base1 + [2:nang 1];

    nloc = numel(i0);
    rows = frow:(frow+nloc-1);
    faces(rows,:) = [i0(:) i1n(:) i1(:)];
    frow = frow + nloc;

    rows = frow:(frow+nloc-1);
    faces(rows,:) = [i0(:) i0n(:) i1n(:)];
    frow = frow + nloc;
end

if capEnds
    % Hemispheres built with their rim on z = 0: the bottom one (pole at
    % z = -r) is already in place at p1, the top one (pole at z = +r) is
    % moved to the end of the axis, z = L.
    [vCapBottom, fCapBottom] = hemisphereMesh(r, nang, nCapRings, -1);
    [vCapTop, fCapTop] = hemisphereMesh(r, nang, nCapRings, 1);

    vCapTop = vCapTop + [0 0 L];

    % The rim of each hemisphere (its last nang vertices) coincides with
    % the first (bottom) or last (top) ring of the side surface, with the
    % same angular ordering. Only the pole and the inner rings are added,
    % and the cap faces are re-indexed so that rim vertices point to the
    % side rings: map(k) is the global index of local hemisphere vertex k.
    nSideVerts = size(vertices,1);
    nCapVerts = 1 + (nCapRings-1)*nang;       % pole + inner rings
    mapBottom = [nSideVerts + (1:nCapVerts), 1:nang];
    mapTop    = [nSideVerts + nCapVerts + (1:nCapVerts), (nz-1)*nang + (1:nang)];

    vertices = [vertices; vCapBottom(1:nCapVerts,:); vCapTop(1:nCapVerts,:)];

    nfBottom = size(fCapBottom,1);
    nfTop = size(fCapTop,1);

    faces(frow:(frow+nfBottom-1),:) = mapBottom(fCapBottom);
    frow = frow + nfBottom;

    faces(frow:(frow+nfTop-1),:) = mapTop(fCapTop);
end

% Rotate local +z axis to (p2-p1) and translate to p1.
% Rodrigues' formula: the rotation by angle vartheta about the unit axis k
% is R = I + sin(vartheta)*K + (1 - cos(vartheta))*K^2, with K the
% cross-product matrix of k (K*x = cross(k, x)). Taking
% k = cross(e_z, v)/|cross(e_z, v)| gives cos(vartheta) = dot(e_z, v) = c
% and sin(vartheta) = |cross(e_z, v)| = s. When v is (anti)parallel to e_z,
% k is undefined and the two cases are handled separately.
v = axisVec / L;
z0 = [0 0 1];
c = dot(z0,v);

if c > 1 - 1e-12
    R = eye(3);
elseif c < -1 + 1e-12
    % 180 deg rotation around x.
    R = [1 0 0; 0 -1 0; 0 0 -1];
else
    k = cross(z0,v);
    s = norm(k);
    k = k / s;

    K = [
         0   -k(3)  k(2);
        k(3)   0   -k(1);
       -k(2)  k(1)   0
    ];

    R = eye(3) + s*K + (1-c)*(K*K);
end

% Apply R to every vertex (rows of vertices), then translate.
vertices = (R * vertices.').';
vertices = vertices + p1;

end

function [vertices, faces] = hemisphereMesh(radius, nang, nRings, signZ)

% Local hemisphere with pole at z = signZ * radius and rim on z = 0.
% Vertex 1 is the pole; ring k (k = 1..nRings, at polar angle
% k*pi/(2*nRings)) occupies rows 1 + (k-1)*nang + (1:nang). Ring nRings is
% the rim. Faces: a fan of nang triangles from the pole to ring 1, then
% 2*nang triangles between each pair of consecutive rings.
phi = (0:nang-1)' * (2*pi/nang);
theta = linspace(0, pi/2, nRings+1).';

vertices = zeros(1 + nRings*nang, 3);
vertices(1,:) = [0 0 signZ * radius];

% Rings at polar angles theta(2:end) (theta(1) = 0 is the pole): ring
% radius r*sin(t), height signZ*r*cos(t).
idx = 2;   % next free row of vertices
for ir = 2:numel(theta)
    t = theta(ir);
    rr = radius * sin(t);
    zz = signZ * radius * cos(t);

    ring = [
        rr * cos(phi), ...
        rr * sin(phi), ...
        zz * ones(nang,1)
    ];
    vertices(idx:idx+nang-1,:) = ring;
    idx = idx + nang;
end

faces = zeros(2*nang*(nRings-1) + nang, 3);
frow = 1;

% Faces are first ordered for the top hemisphere (signZ > 0), where they
% are anticlockwise seen from outside (outward normals).
% Fan from pole to first ring.
ring1 = 2:(1+nang);
ring1n = [3:(1+nang) 2];
faces(frow:frow+nang-1,:) = [ones(nang,1) ring1(:) ring1n(:)];
frow = frow + nang;

% Band between rings ir-1 and ir (ring k starts at row 2 + (k-1)*nang):
% two triangles per angular segment, as for the cylinder sides.
for ir = 2:nRings
    base0 = 2 + (ir-2)*nang;   % first vertex of ring ir-1
    base1 = base0 + nang;      % first vertex of ring ir

    i0 = base0 + (0:nang-1);
    i1 = base1 + (0:nang-1);
    i0n = base0 + [1:nang-1 0];
    i1n = base1 + [1:nang-1 0];

    faces(frow:frow+nang-1,:) = [i0(:) i1(:) i1n(:)];
    frow = frow + nang;
    faces(frow:frow+nang-1,:) = [i0(:) i1n(:) i0n(:)];
    frow = frow + nang;
end

% The bottom hemisphere is the mirror image (z -> -z) of the top one, and a
% reflection reverses orientation: swap two vertices of every face to keep
% the normals outward.
if signZ < 0
    faces = faces(:, [1 3 2]);
end

end
