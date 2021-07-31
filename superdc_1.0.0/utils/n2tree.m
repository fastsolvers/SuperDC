function tr = n2tree(n)
% n: number of nodes of a full binary tree
% The number of leaves is (n+1)/2
% Result: for any odd n, there always exists a full binary tr

if mod(n,2) == 0
    error('input must be odd');
end

if n == 1
    tr = 0;
    return;
else
    n1 = 2^floor(log2(n))-1;
    tr1 = btree(n1);
    tr2 = n2tree(n-n1-1)+n1;
    tr1(end) = n; tr2(end) = n;
    tr = [tr1 tr2 0];
end