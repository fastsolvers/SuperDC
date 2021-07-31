function [p,l] = btree(n)
% Binary tree.
% p(j) is the parent of column j in the tree, or 0 if j is a root.

if mod(n,2) == 0
    error('input must be odd');
end

if n > 3
    m = (n-1)/2;
    [t1,l1] = btree(m);
    t2 = m+t1;
    t1(end) = n;
    t2(end) = n;
    p = [t1 t2 0];
    l = [l1+1 l1+1 1];
elseif n == 3
    p = [3 3 0];
    l = [2 2 1];
else
    p = 0;
    l = 1;
end