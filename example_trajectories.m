% Example wall
clearvars
addpath(genpath('./'))

params_reference

params.N_cells = 50;
params.deltat_save = params.deltat;
wall.amplitudeOsc = 0.5;
angle = 0;

[x_pts, y_pts, fibre_cols] = wall_base(cfg, wall, angle);
[F, V] = wall_to_triangles(x_pts, y_pts, fibre_cols, cfg, wall);
[centres, normals] = triangles_to_centres(F,V);
[X, t_exit] = movement_scaffold(params, limits, centres, normals, initials);

plot_trajectories(X, limits, params, F, V)

save("results\wall_trajectory", "X", "limits", "params", "F", "V")

%%
% Example fibres

clearvars
addpath(genpath('./'))

params_reference

params.N_cells = 50;
params.deltat_save = params.deltat;
wall.amplitudeOsc = 0.5;
angle = 0;

mesh =  mesh_fibres_aligned(cfg, angle);
[F,V] = mesh_to_triangles(mesh, fib);
[centres, normals] = triangles_to_centres(F,V);
[X, t_exit] = movement_scaffold(params, limits, centres, normals, initials);

plot_trajectories(X, limits, params, F, V)

save("results\fibres_trajectory", "X", "limits", "params", "F", "V")

%%
% Example Voronoi with low coefficient
clearvars
addpath(genpath('./'))

params_reference

params.N_cells = 50;
params.deltat_save = params.deltat;

rng(2)
Npoints = 15;
a = 0.5;

points = [
    cfg.xmin + (cfg.xmax - cfg.xmin) * rand(Npoints, 1), ...
    cfg.ymin + (cfg.ymax - cfg.ymin) * rand(Npoints, 1), ...
    cfg.zmin + (cfg.zmax- cfg.zmin) * rand(Npoints, 1)
];


A = [1 0 0; 0 a 0; 0 0 a];
points_edited = points * A;
cfg_tmp = cfg;
tmp = num2cell([cfg.xmax cfg.ymax cfg.zmax] * A);
[cfg_tmp.xmax, cfg_tmp.ymax, cfg_tmp.zmax] = tmp{:};
tmp = num2cell([cfg.xmin cfg.ymin cfg.zmin] * A);
[cfg_tmp.xmin, cfg_tmp.ymin, cfg_tmp.zmin] = tmp{:};

mesh = mesh_voronoi(points_edited, cfg_tmp);
    
mesh.vertices = mesh.vertices / A;
mesh.points = mesh.points / A;

[F, V] = mesh_to_triangles(mesh, fib);

[centres, normals] = triangles_to_centres(F,V);
[X, t_exit] = movement_scaffold(params, limits, centres, normals, initials);

plot_trajectories(X, limits, params, F, V)

save("results\voronoi_low_trajectory", "X", "limits", "params", "F", "V", "a")

% TR = triangulation(F, V);
% stlwrite(TR, stlFile);
