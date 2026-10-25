function plot_trajectories(X, limits, params, F, V)
% PLOT_TRAJECTORIES Plots the trajectories of cells in 3D space along with
% the scaffold if provided.

% figure

if ~isempty(X)
    N_cells = size(X,1);

    for k=1:N_cells
        xx=X(k,1,:);
        yy=X(k,2,:);
        zz=X(k,3,:);
        plot3(xx(:),yy(:),zz(:),'LineWidth',2)
        % plot3(X(1,k,:),X(2,k,:),X(3,k,:),'LineWidth',2)
        hold on
    end

    [e,f,g]=sphere;
    e1=e*params.r+X(1,1,1);
    f1=f*params.r+X(1,2,1);
    g1=g*params.r+X(1,3,1);
    surf(e1,f1,g1,'EdgeColor', 'none', 'FaceColor', 0.5*[1 1 1])

end


if nargin > 3 && ~isempty(F) && ~isempty(V)
    trisurf(F, V(:,1), V(:,2), V(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'none');
end

if nargin < 4 && ~strcmp(params.scaffold,'no_scaffold') && ~strcmp(params.scaffold,' ')
    [stlStruct,~] = stlread("scaffolds/" + params.scaffold + ".stl");
    V = stlStruct.Points;
    F = stlStruct.ConnectivityList;
    
    trisurf(F, V(:,1), V(:,2), V(:,3), 'FaceColor', 'cyan', 'EdgeColor', 'none');
end

camlight headlight;
lighting gouraud;
axis equal
xlim([limits.xmin limits.xmax])
ylim([limits.ymin limits.ymax])
zlim([limits.zmin limits.zmax])