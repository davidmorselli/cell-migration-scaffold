% Example wall
clearvars
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
%%
% Example other scaffold
clearvars
addpath(genpath('./'))

params_reference
params.scaffold = 'SSW_random_Radius05';

% params.rot_exp = 2; % Exponent of rotation matrix for the scaffold. The default is no rotation
% params.rot_mat = [0 -1 0; 1 0 0; 0 0 1]^params.rot_exp; % Actual rotation matrix


[stlStruct,stlName] = stlread("scaffolds/"+params.scaffold + ".stl");
V = stlStruct.Points;
F = stlStruct.ConnectivityList;
V = V * params.rot_mat';

limits={};
limits.xmax = max(V(:,1)) + params.xmaxeff;
limits.xmin = min(V(:,1)) - params.xmineff;
limits.ymax = max(V(:,2));
limits.ymin = min(V(:,2));
limits.zmax = max(V(:,3));
limits.zmin = min(V(:,3));

initials = [limits.xmin + params.r 
    (limits.ymax + limits.ymin) / 2 
    (limits.zmax + limits.zmin) / 2 ];

params.N_cells = 10;
params.T_macro = 50;
params.deltat_save = params.deltat*10;

[centres, normals] = triangles_to_centres(F,V);
tic
[X, t_exit] = movement_scaffold(params, limits, centres, normals, initials);
toc

figure
plot_trajectories(X, limits, params, F, V)

save("results\"+string(params.scaffold)+"trajectories", "X", "limits", "params", "F", "V")
%%
% Example other scaffold
clearvars
addpath(genpath('./'))

params_reference
params.scaffold = 'shg';

params.rot_exp = 2; % Exponent of rotation matrix for the scaffold. The default is no rotation
params.rot_mat = [0 -1 0; 1 0 0; 0 0 1]^params.rot_exp; % Actual rotation matrix


[stlStruct,stlName] = stlread("scaffolds/"+params.scaffold + ".stl");
V = stlStruct.Points;
F = stlStruct.ConnectivityList;
V = V * params.rot_mat';

limits={};
limits.xmax = max(V(:,1)) + params.xmaxeff;
limits.xmin = min(V(:,1)) - params.xmineff;
limits.ymax = max(V(:,2));
limits.ymin = min(V(:,2));
limits.zmax = max(V(:,3));
limits.zmin = min(V(:,3));

initials = [-120 %limits.xmin + params.r 
    (limits.ymax + limits.ymin) / 2 
    (limits.zmax + limits.zmin) / 2 ];

params.N_cells = 10;
params.T_macro = 200;
params.b=0.4;
params.deltat_save = params.deltat*10;

[centres, normals] = triangles_to_centres(F,V);
tic
[X, t_exit] = movement_scaffold(params, limits, centres, normals, initials);
toc

figure
plot_trajectories(X, limits, params, F, V)

save("results\"+string(params.scaffold)+"trajectories", "X", "limits", "params", "F", "V")
%%

clearvars
params_reference

params.N_cells = 50;
params.deltat_save = params.deltat;

rng(2)
Npoints = 15;
a = 2;

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

figure
plot_trajectories(X, limits, params, F, V)

save("results\voronoi_high_trajectory", "X", "limits", "params", "F", "V", "a")




