function d = crossdiff(a,ind)
%Function that computes the cross product of the difference of vectors and
%then normalises it; a has size n times 3 and ind is a vector (only the
%first 4 elements are used) that indicates which rows of a should be used.

d=cross(a(ind(1),:)-a(ind(3),:),a(ind(2),:)-a(ind(4),:));
n=norm(d);

if n~=0
    d=d/norm(d);
end
end