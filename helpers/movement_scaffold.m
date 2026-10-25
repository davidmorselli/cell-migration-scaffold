function [X, t_exit] = movement_scaffold(params, limits, centres, normals, initial)
%Contact inhibition for cells closer than rmin
%   MOST ADVANCED VERSION SO FAR AND WRITTEN IN A CLEAR WAY



T = params.T_macro / params.deltat;
v = params.v_macro * params.deltat;
prob_change_free = params.rate_change_free * params.deltat;
prob_change_bound = params.rate_change_bound * params.deltat;
prob_bind = params.rate_bind * params.deltat;
prob_unbind = params.rate_unbind * params.deltat;

save_stride = params.deltat_save / params.deltat;

X = nan(params.N_cells, 7, params.T_macro / params.deltat_save + 1); %x_1, x_2, x_3, v_1, v_2, v_3, bind

t_exit = nan(params.N_cells,1);


% Set random seed
if ~isnan(params.rand)
    rng(params.rand, 'twister') %the generator type is added to ensure reproducibility inside parfor
end

tree = KDTreeSearcher(centres);

% Initialise positions
if nargin >= 5 && ~isempty(initial)
    X(:,1,1) = initial(1);
    X(:,2,1) = initial(2);
    X(:,3,1) = initial(3);
else
    X(:,1,1) = limits.xmin;
    % X(:,1,1) = (limits.xmin+limits.xmax) / 2 - (- limits.xmin+ limits.xmax)/4;
    X(:,2,1) = (limits.ymax + limits.ymin) / 2;
    X(:,3,1) = (limits.zmax + limits.zmin) / 2;
end
X(:,4,1) = v;
X(:,5:7,1) = 0;

% Ensure binding state at t=0 is zero (unbound)
X(:,7,1) = 0;

% Initialize status logical arrays
still = false(params.N_cells, 1);
new_exit = false(params.N_cells, 1);
Y_prec = X(~still, :, 1);

% Find closest scaffold point and normal
[ind, dist] = knnsearch(tree, Y_prec(:, 1:3));
near = dist <= params.r;
a = normals(ind, :);
a(~near, :) = NaN;
    

% Simulate cell movement
for i=1:T

    if all(still)
        break
    end


    % Select cells that keep moving
    Y = nan(sum(~still), 7);


    Y(~near, 7) = 0; % Cells far from collagen cannot be bound

    % Previously bounded
    % If bound, then randomly unbound
    potential_unbound = (Y_prec(:, 7)==1) & near;
    Y(potential_unbound, 7) = ~(rand(sum(potential_unbound), 1) < prob_unbind); % It becomes 0 with probability prob_unbind


    % Previously unbounded
    % If collagen is close, then randomly bounds
    potential_bound = (Y_prec(:, 7) == 0) & near;
    Y(potential_bound, 7) = rand(sum(potential_bound), 1) < prob_bind;

    now_bound = Y(:,7);


    % Movement

    % Now bounded, change direction
    change = now_bound & (rand(size(Y, 1), 1) < prob_change_bound);
    n_change = sum(change);
    [xx,yy,zz] = sph2cart(2 * pi * rand(n_change, 1), asin(2 * rand(n_change, 1) - 1), 1);
    newdir = (1-params.b) * [xx, yy, zz] + params.b * params.d_bias;
    Y(change, 4:6)  = newdir - sum(newdir .* a(change, :), 2) .* a(change, :); % Project direction tangential to scaffold 
    % Now bounded, same direction reprojected
    newdir = Y_prec((now_bound & ~change), 4:6);
    Y((now_bound & ~change), 4:6) = newdir - sum(newdir .* a(now_bound & ~change, :), 2) .* a(now_bound & ~change, :);


    % Now unbounded, change direction
    change = ~now_bound & (rand(size(Y, 1), 1) < prob_change_free);
    n_change = sum(change);
    [xx,yy,zz] = sph2cart(2 * pi * rand(n_change, 1),asin(2 * rand(n_change, 1) - 1),1);
    Y(change, 4:6) = (1 - params.b)*[xx, yy, zz] + params.b * params.d_bias;
    % Now unbounded, same direction
    Y((~now_bound & ~change), 4:6) = Y_prec((~now_bound & ~change), 4:6);


    % Normalise velocity and update positions
    Y(:, 4:6) = Y(:, 4:6) ./ vecnorm(Y(:, 4:6), 2, 2) * v;
    Y(:, 1:3) = Y_prec(:, 1:3) + Y(:, 4:6);

    % Bouncing from boundaries (this does not take into account chemotaxis)
    going_out = Y(:, 1) < limits.xmin | ...
        Y(:, 2) < limits.ymin | Y(:, 2) > limits.ymax |... 
        Y(:, 3) < limits.zmin | Y(:, 3) > limits.zmax; 

    if any(going_out)
        Y(going_out, 1:6) = BoundaryReflection(Y(going_out, 1:6), Y_prec(going_out, 1:6), limits);
    end

    % Bouncing from fibres
    [ind, dist] = knnsearch(tree, Y(:,1:3));
    near = dist <= params.r;
    a = normals(ind, :);
    a(~near, :) = NaN;
    
    hitting = dist <= params.rmin;

    if any(hitting)
        [Y(hitting, 1:6), near(hitting), a(hitting,:)] = updateHitting( ...
            Y(hitting, 1:6), Y_prec(hitting, 1:6), dist(hitting), ind(hitting), ...
            near(hitting), a(hitting,:), normals, centres, limits, params.r, tree);
    end    


    % Save positions
    if mod(i, save_stride) == 0
       X(~still, :, i / save_stride + 1) = Y;
    end

    arrived = Y(:,1) > limits.xmax; 
    new_exit(still) = false;
    new_exit(~still) = arrived;
    t_exit(new_exit) = i * params.deltat;
    
    still(~still) = arrived; % Take into account cells that need to stop now

    % Prepare for the next iteration
    near = near(~arrived);
    a = a(~arrived,:);
    Y_prec = Y(~arrived, :);

end
end

%------------------------------------------------------------------

function Y = BoundaryReflection(Y, Y_prec, limits)

    xHit = Y(:,1) < limits.xmin;

    yHit = Y(:,2) < limits.ymin | ...
           Y(:,2) > limits.ymax;

    zHit = Y(:,3) < limits.zmin | ...
           Y(:,3) > limits.zmax;

    Y(xHit,4) = -Y(xHit,4);
    Y(yHit,5) = -Y(yHit,5);
    Y(zHit,6) = -Y(zHit,6);

    Y(:,1:3) = Y_prec(:,1:3) + Y(:,4:6);

end

function [y, v, ind, dist] = HitAndOut(temp, temp_prec, y, ...
                        c_min, a, limits, tree)

    N = size(y,1);

    % Boundary normals
    n = zeros(N,3);

    n(:,1) = n(:,1) + (y(:,1) < limits.xmin);

    n(:,2) = n(:,2) + (y(:,2) < limits.ymin);
    n(:,2) = n(:,2) - (y(:,2) > limits.ymax);

    n(:,3) = n(:,3) + (y(:,3) < limits.zmin);
    n(:,3) = n(:,3) - (y(:,3) > limits.zmax);

    % Normal from the scaffold
    d = temp(:,1:3) - c_min; %temp is the point that was rejected because too close to the scaffold
    a = a .* sign(sum(d .* a, 2)); %this ensures that a points outwards

    % Combined reflection direction
    vect = a ./ vecnorm(a,2,2) + n ./ vecnorm(n,2,2);

    vect = vect ./ vecnorm(vect,2,2);

    % Reflect velocity
    v = temp(:,4:6) - 2 * sum(temp(:,4:6) .* vect, 2) .* vect;

    % Update positions
    y = temp_prec(:,1:3) + v;

    % Nearest-neighbour search
    [ind, dist] = knnsearch(tree, y);

end


%------------------------------------------------------------------

function [Y, near, a] = updateHitting(Y, Y_prec, dist, indtemp, near, a, normals, centres, limits, r, tree)
%NB: Y = Y(hitting, 1:6)

% Reflection step
v = Y(:, 4:6) - 2 * sum(Y(:, 4:6) .* a, 2) .* a;
y = Y_prec(:, 1:3) + v;

% Nearest scaffold check
[ind_new, dist_new] = knnsearch(tree, y);
getting_closer = dist_new < dist;

% If reflection brings cell closer, then it
% should not reflect and the speed should be the one previously computed
v(getting_closer, :) = Y(getting_closer, 4:6);
y = Y_prec(:, 1:3) + v;

% Boundary check
hit_and_out = y(:,1) < limits.xmin | ...
              y(:,2) < limits.ymin | ...
              y(:,2) > limits.ymax | ...
              y(:,3) < limits.zmin | ...
              y(:,3) > limits.zmax;

if any(hit_and_out)
    [y(hit_and_out,:), v(hit_and_out,:), ind_new(hit_and_out), dist_new(hit_and_out)] = ...
        HitAndOut(Y(hit_and_out,1:6), ...
                  Y_prec(hit_and_out,1:6), ...
                  y(hit_and_out,:), ...
                  centres(indtemp(hit_and_out),:), ...
                  a(hit_and_out,:), ...
                  limits, tree); % the centres are the ones close to the position Y that was hitting
end

% Update variables
a(~getting_closer, :) = normals(ind_new(~getting_closer), :);
near(~getting_closer) = dist_new(~getting_closer) <= r;
Y(:, 1:3)   = y;
Y(:, 4:6)   = v;

end