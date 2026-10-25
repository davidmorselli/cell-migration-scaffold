function [F, V] = wall_to_triangles(x_pts, y_pts, fibre_cols, cfg, wall)

%WALL_TO_TRIANGLES  Extrude the 2-D wall outlines from wall_base along z and
%                   tessellate their lateral surface into a triangle mesh,
%                   ready for use with MATLAB's triangulation() or stlwrite().
%
%   [F, V] = wall_to_triangles(x_pts, y_pts, fibre_cols, cfg, wall)
%
%   Every outline point is copied onto wall.nZLayer equally spaced layers
%   z = linspace(cfg.zmin, cfg.zmax, wall.nZLayer). Each pair of
%   neighbouring outline points (along the orig line, along the shift line
%   and around the two end-cap arcs) is joined by a strip of 2*(nZ-1)
%   triangles. Only the lateral surface is built: the walls are open at
%   z = cfg.zmin and z = cfg.zmax (no top/bottom faces).
%   The faces are oriented so that their normals point out of each wall:
%   the outline is walked clockwise (seen from +z), i.e. shift line
%   left -> right, right cap shift -> orig, orig line right -> left, left
%   cap orig -> shift (implemented by flipping the orig and right-cap
%   strips, see appendStrip).
%
%   Inputs
%     x_pts      – N×1 outline point x-coordinates (from wall_base)
%     y_pts      – N×1 outline point y-coordinates (from wall_base)
%     fibre_cols – connectivity struct from wall_base
%     cfg        – scaffold struct (from params_reference): zmin, zmax
%     wall       – wall struct (from params_reference): nZLayer
%
%   Outputs
%     F – M×3 face array    (vertex indices, one triangle per row),
%         M = 2*(nZ-1)*N
%     V – K×3 vertex array  (x y z coordinates, one vertex per row),
%         K = N*nZ
%
%   Typical usage (see build_parallel_walls.m):
%     [x_pts, y_pts, fibre_cols] = wall_base(cfg, wall, angle);
%     [F, V] = wall_to_triangles(x_pts, y_pts, fibre_cols, cfg, wall);
%     TR = triangulation(F, V);
%     stlwrite(TR, 'walls.stl');

    nPts  = numel(x_pts);
    z     = linspace(cfg.zmin, cfg.zmax, wall.nZLayer);
    nZ    = wall.nZLayer;
    nCols = numel(fibre_cols(1).straight) / 2;

    % Build vertex array
    % Layout: vertex index = (pt-1)*nZ + zi
    %   pt  → 2-D outline point index (1..nPts)
    %   zi  → z-layer index           (1..nZ)
    %
    %   V = [ x1,y1,z1;  x1,y1,z2;  ...  x1,y1,zEnd;
    %         x2,y2,z1;  x2,y2,z2;  ...  x2,y2,zEnd;
    %         ...                                     ]

    V = [ repelem(x_pts(:), nZ), ...
          repelem(y_pts(:), nZ), ...
          repmat(z(:), nPts, 1) ];

    % Exact number of faces: each outline is a closed loop, so there is one
    % strip (2*(nZ-1) triangles) per outline point.
    nFacesMax = 2 * (nZ - 1) * nPts;
    F  = zeros(nFacesMax, 3);
    nF = 0;   % faces stored so far (shared with the nested functions)

    % Index in V of outline point pt on z-layer zi (see layout above).
    vid = @(pt, zi) (pt - 1) * nZ + zi;

    % Store one triangle in the next free row of F.
    function appendFace(a, b, c)
        nF = nF + 1;
        F(nF,:) = [a, b, c];
    end

    % Join outline points ptA and ptB across all z-layers with 2*(nZ-1)
    % triangles; flip reverses the vertex order (i.e. the normal direction).
    % Normal for flip = false: e_z x (P_B - P_A), i.e. P_B - P_A rotated by
    % +90 degrees (to the left of the direction of travel); for
    % flip = true it points to the right.
    function appendStrip(ptA, ptB, flip)
        % Each quad between layers i and i+1 is split into two triangles.
        % Naming: v<layer offset><point>, e.g. v10 = point A on layer i+1.
        for i = 1:nZ-1
            v00 = vid(ptA, i);
            v10 = vid(ptA, i+1);
            v01 = vid(ptB, i);
            v11 = vid(ptB, i+1);

            if ~flip
                appendFace(v00, v10, v01);
                appendFace(v11, v01, v10);
            else
                appendFace(v00, v01, v10);
                appendFace(v11, v10, v01);
            end
        end
    end

    % One closed loop of strips per wall.
    for r = 1:numel(fibre_cols)
        st = fibre_cols(r).straight;
        orig = st(1:nCols);
        shift = st(nCols+1:2*nCols);

        % Orig line (lower edge): travelling in +x, the outside (-y) is to
        % the right, so flip.
        for c = 1:nCols-1
            appendStrip(orig(c), orig(c+1), true);
        end

        % Shift line (upper edge): travelling in +x, the outside (+y) is to
        % the left, so no flip.
        for c = 1:nCols-1
            appendStrip(shift(c), shift(c+1), false);
        end

        % Left/right end caps close off the stadium-shaped cross-section;
        % the two long edges (orig/shift) are only joined here, at the
        % ends — no cross-links should exist along the wall's length.
        % arc_left is stored ordered from near-shift to near-orig, so
        % reverse it to walk orig(1) -> ... -> shift(1) via the closest
        % neighbouring arc points instead of jumping across the cap.
        arc_left_ordered = fliplr(fibre_cols(r).arc_left(:)');
        % Walking orig(1) -> shift(1) around the left cap, the outside (-x)
        % is to the left: no flip.
        strip_L = [orig(1), arc_left_ordered, shift(1)];
        for j = 1:numel(strip_L)-1
            appendStrip(strip_L(j), strip_L(j+1), false);
        end

        % Walking orig(end) -> shift(end) around the right cap, the outside
        % (+x) is to the right: flip.
        strip_R = [orig(end), fibre_cols(r).arc_right(:)', shift(end)];
        for j = 1:numel(strip_R)-1
            appendStrip(strip_R(j), strip_R(j+1), true);
        end
    end

end