% Parameters for simulations

params.rand = 1; % Random seed
params.N_cells = 1000; % Number of cells
params.T_macro = 1000; % Total time [min]
params.deltat = 0.01; % Time discretisation [min]
params.deltat_save = 2; % Time discretisation [min]

params.v_macro = 5; % Speed [micron/min]
params.r = 3.5; % Interaction radius of cells with fibres [micron]
params.rmin = 2; % Minimal distance between cells and fibres [micron]
params.rate_change_free = 0.75; % Rate of direction change when unbounded [1/min]
params.rate_change_bound = 0.75; % Rate of direction change when bounded [1/min]
params.rate_bind = 2;  % Binding rate of cells to scaffold [1/min]
params.rate_unbind = 0.6; % Unbinding rate of cells from scaffold [1/min]

params.b = 0.3; % Bias strength
params.d_bias = [1 0 0]; % Bias direction
params.xmaxeff = 0; % Additional distance from maximal x cells are allowed to go
params.xmineff = 20; % Additional distance from minimal x cells are allowed to go

params.rot_exp = 0; % Exponent of rotation matrix for the scaffold. The default is no rotation
params.rot_mat = [0 -1 0; 1 0 0; 0 0 1]^params.rot_exp; % Actual rotation matrix

params.scaffold = ' '; % The name of the scaffold is only present if the scaffold has been saved

% Limits for the simulation domain

limits.xmax = 50;
limits.xmin = -20;
limits.ymax = 50;
limits.ymin = 0;
limits.zmax = 30;
limits.zmin = 0;

% Initial condition

initials = [limits.xmin + params.r 
    (limits.ymax + limits.ymin) / 2 
    (limits.zmax + limits.zmin) / 2 ];

% Parameters for scaffolds
% Fibres
fib.cylRadius = 0.5;
fib.nang = 16;
fib.nCapRings = 5;
fib.ringSpacing = 0.5;
fib.capEnds = true;

% Walls
wall.size = 2 * fib.cylRadius;
wall.nCapRings = 2 * fib.nCapRings +1;
wall.dx = fib.ringSpacing;
wall.nZLayer = 30;
wall.amplitudeOsc = 0;
wall.periodOsc = 10;
wall.phaseOsc = 0;

% Scaffold dimension
cfg.xmin = 0 + fib.cylRadius;
cfg.xmax = 50 - fib.cylRadius; 
cfg.ymin = 0 + fib.cylRadius;
cfg.ymax = 50 - fib.cylRadius;
cfg.zmin = 0 + fib.cylRadius;
cfg.zmax = 30 - fib.cylRadius;

% cfg.dx = 0.5; %ANALOGO TO RINGSPACING, CONTROLLARE

cfg.yspace = 12; % only relevant for parallel walls and fibres
cfg.zspace = 9.5; % only relevant for parallel fibres
