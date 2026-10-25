function [centres, normals] = triangles_to_centres(F,V)
% Compute baricentres and normals from triangulation

centres = zeros(size(F));
normals = centres;

for i=1:size(F,1)
    centres(i,:)=(V(F(i,1),:)+V(F(i,2),:)+V(F(i,3),:))/3;
    normals(i,:)=crossdiff(V,[F(i,:) F(i,3)]);
end