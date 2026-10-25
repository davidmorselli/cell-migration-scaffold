function [x_pts, y_pts, fibre_cols] = wall_base(cfg, wall, alpha)
%WALL_BASE  2-D cross-section (in the x-y plane) of a set of parallel wavy
%           walls, to be extruded in z by wall_to_triangles.
%   [x_pts, y_pts, fibre_cols] = wall_base(cfg, wall, alpha)
%
%   Each wall is a strip of thickness wall.size in y whose lower edge
%   ("orig" line) follows
%     y = y0 + f(x),  f(x) = wall.amplitudeOsc*sin(2*pi*(x - cfg.xmin)/wall.periodOsc + wall.phaseOsc),
%   sampled at x = cfg.xmin:wall.dx:cfg.xmax, with
%     y0 = (cfg.ymin - wall.size/2 + wall.amplitudeOsc) : cfg.yspace : cfg.ymax,
%   so that the lowest point of the first wall is at cfg.ymin - wall.size/2
%   (y = 0 with the values in params_reference).
%   The upper edge ("shift" line) is the lower one translated by wall.size
%   in y. The two edges are joined at each end by a semicircular arc of
%   radius wall.size/2 (stadium-shaped outline).
%   The walls are built at alpha = 0 and then each wall is rotated rigidly
%   by alpha about its own left end-cap centre (as the fibres of
%   mesh_fibres_aligned rotate about the centre of their left cap), so for
%   alpha ~= 0 the walls no longer end at cfg.xmax and may extend beyond
%   cfg.ymax.
%
%   Inputs
%     cfg   - scaffold struct (from params_reference): xmin, xmax, ymin,
%             ymax, yspace
%     wall  - wall struct (from params_reference): size (thickness),
%             dx (sampling step in x), nCapRings (number of points on each
%             end-cap semicircle, including its two end points, which are
%             the end points of the orig/shift lines), amplitudeOsc,
%             periodOsc, phaseOsc
%     alpha - rotation angle [degrees], anticlockwise
%
%   Outputs
%     x_pts, y_pts - Nx1 coordinates of all outline points of all walls
%     fibre_cols   - struct array (one entry per wall) with the indices
%                    into x_pts/y_pts of
%                      straight  : orig line (columns 1..nCols) followed by
%                                  shift line (columns 1..nCols)
%                      arc_left  : interior points of the left cap, ordered
%                                  from the shift side to the orig side
%                      arc_right : interior points of the right cap, ordered
%                                  from the orig side to the shift side
%                    Per wall, the points are stored in this same order
%                    (2*nCols + 2*(wall.nCapRings - 2) points).
%
%   Typical usage (see build_parallel_walls.m):
%     [x_pts, y_pts, fibre_cols] = wall_base(cfg, wall, angle);
%     [F, V] = wall_to_triangles(x_pts, y_pts, fibre_cols, cfg, wall);


    % Wall profile (sinusoidal oscillation in y)
    f = @ (x) wall.amplitudeOsc*sin(2*pi*(x-cfg.xmin)/wall.periodOsc + wall.phaseOsc);

    % Build alpha = 0 geometry in local 2-D
    % Common x-sampling and profile of all walls; wall r is the profile
    % shifted by yStarts(r).
    x_base = cfg.xmin:wall.dx:cfg.xmax;
    y_base = f(x_base);

    % yStarts = (cfg.ymin - wall.size / 2 + wall.amplitudeOsc) : cfg.yspace : cfg.ymax;
    yStarts = (cfg.ymin - wall.size / 2) : cfg.yspace : cfg.ymax;

    % Sizes: nCols points on each straight line, nInterior points inside
    % each end-cap arc (its two end points are the line end points).
    nRows         = numel(yStarts);
    nCols         = numel(x_base);
    nInterior     = wall.nCapRings - 2;
    nPerRow       = 2*nCols + 2*nInterior;
    nTotal        = nRows * nPerRow;

    x_pts   = zeros(nTotal, 1);
    y_pts   = zeros(nTotal, 1);

    fibre_cols = struct('straight', cell(nRows,1), ...
                        'arc_left', cell(nRows,1), ...
                        'arc_right', cell(nRows,1));

    % Rotate each wall by alpha about its own left end-cap centre
    % R is the anticlockwise rotation by alpha in the x-y plane.
    theta = deg2rad(alpha);
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

    idx = 0;   % number of points stored so far
    for r = 1:nRows
        y_orig_row  = y_base + yStarts(r);
        y_shift_row = y_orig_row + wall.size;

        % Straight segments in alpha=0 frame.
        x_row = [x_base, x_base].';
        y_row = [y_orig_row, y_shift_row].';

        % End-cap arcs in alpha=0 frame (interior points only).
        % Left arc goes from the shift side (theta = pi/2) to the orig side
        % (theta = 3*pi/2); right arc from orig (-pi/2) to shift (pi/2).
        % Each arc is a semicircle centred midway between the orig and shift
        % end points, with radius half their distance (= wall.size/2).
        yc_left = 0.5 * (y_orig_row(1) + y_shift_row(1));
        rArc    = 0.5 * (y_shift_row(1) - y_orig_row(1));
        yc_right = 0.5 * (y_shift_row(end) + y_orig_row(end));

        % wall.nCapRings equally spaced angles; the first and last are the
        % line end points and are dropped.
        thetaL = linspace( pi/2, 3*pi/2, wall.nCapRings);
        thetaR = linspace(-pi/2,   pi/2, wall.nCapRings);
        tL = thetaL(2:end-1);
        tR = thetaR(2:end-1);

        x_arc_left  = x_base(1)   + rArc * cos(tL);
        y_arc_left  = yc_left     + rArc * sin(tL);

        x_arc_right = x_base(end) + rArc * cos(tR);
        y_arc_right = yc_right    + rArc * sin(tR);

        % All points of this wall, in the order orig, shift, left arc,
        % right arc.
        x0 = [x_row; x_arc_left(:); x_arc_right(:)];
        y0 = [y_row; y_arc_left(:);  y_arc_right(:)];

        % Rotation about the centre of the left end cap: translate the pivot
        % to the origin, rotate (row vectors, so multiply by R.'), translate
        % back.
        pivot = [x_base(1), yc_left];

        pts = [x0 - pivot(1), y0 - pivot(2)] * R.';

        rows = idx + (1:numel(x0));
        x_pts(rows)   = pts(:,1) + pivot(1);
        y_pts(rows)   = pts(:,2) + pivot(2);

        % Connectivity for wall_to_triangles: straight block is the
        % first 2*nCols entries (orig then shift, each contiguous);
        % the remaining entries split evenly into left/right cap arcs.
        fibre_cols(r).straight  = rows(1 : 2*nCols);
        fibre_cols(r).arc_left  = rows(2*nCols + 1 : 2*nCols + nInterior);
        fibre_cols(r).arc_right = rows(2*nCols + nInterior + 1 : 2*nCols + 2*nInterior);

        idx = idx + numel(x0);
    end

end

