function mesh = mesh_fibres_aligned(cfg, alpha)
%MESH_FIBRES_ALIGNED Edge mesh (centrelines) of a grid of parallel fibres.
%   mesh = mesh_fibres_aligned(cfg, alpha)
%
%   Builds one straight fibre for every pair (y0, z0), with
%   y0 = cfg.ymin:cfg.yspace:cfg.ymax and z0 = cfg.zmin:cfg.zspace:cfg.zmax.
%   Each fibre starts at (cfg.xmin, y0, z0) and ends at
%   (cfg.xmin + (cfg.xmax - cfg.xmin) * cosd(alpha), 
%   y0 + (cfg.xmax - cfg.xmin)*sind(alpha), z0), 
%   i.e. fibres are rotated.
%
%   Inputs
%     cfg   - scaffold struct (from params_reference) with fields
%             xmin, xmax, ymin, ymax, zmin, zmax (centreline limits),
%             yspace, zspace (spacing between fibres in y and z)
%     alpha - inclination parameter [degrees]
%
%   Output
%     mesh.vertices - (2*nFibres)x3 fibre end points; rows 2k-1 and 2k are
%                     the start and end of fibre k (y varies fastest, then z)
%     mesh.edges    - nFibres x 2 vertex indices [start end] of each fibre
%
%   Typical usage (see build_fibres_aligned.m):
%     mesh  = mesh_fibres_aligned(cfg, angle);
%     [F,V] = mesh_to_triangles(mesh, fib);

% Generate fibre start positions in the y-direction
yStarts = cfg.ymin:cfg.yspace:cfg.ymax;
% (x, y) of the end points of the fibres in one z-plane: rows 2k-1 and 2k
% are the start (cfg.xmin, y0) and end of fibre k. The end is the start
% plus L*(cosd(alpha), sind(alpha)), with L = cfg.xmax - cfg.xmin, i.e.
% each fibre is rotated by alpha about its start point.
xy_matrix = [repmat([cfg.xmin; cfg.xmin + (cfg.xmax - cfg.xmin) * cosd(alpha)], ...
    length(yStarts),1) repelem(yStarts',2)];
xy_matrix(2:2:end,2) = xy_matrix(2:2:end,2) + (cfg.xmax - cfg.xmin) * sind(alpha);
zStarts = (cfg.zmin : cfg.zspace : cfg.zmax)';

% Build mesh: repeat the (x, y) layout for each z-plane; each fibre's edge
% joins its consecutive start/end vertices (1-2, 3-4, ...).
mesh.vertices = [repmat(xy_matrix, length(zStarts), 1) repelem(zStarts, 2*length(yStarts))];
start_vert = 1:2:size(mesh.vertices,1);
mesh.edges = [start_vert' start_vert' + 1];

end 
