%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%  convert a banded matrix into block-tri-diagonal form  %%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [D, F, G] = banded2blktri(A, r)
%%%% Input:
%%%% convert a banded matrix A to a block-tridiagonal representation w/
%%%% block size r, which must be greater than the bandwidth of A

%%% Output:
%%% D: cell array of diagonal blocks
%%% F: cell array of (-1) subdiagonal blocks
%%% G: cell array of (+1) superdiagonal blocks

if nargin < 2
    r = min(2^9, floor(size(A,1)/16));
end

[~, d] = spdiags(A);
bw = max(abs(d));     % bandwidth of A
if r < bw
    r = bw + 1;
end

n = size(A, 1);
[~, m] = npart(n, r);
k = length(m);

D = cell(k, 1);
F = cell(k, 1);
G = cell(k, 1);

l = 1;
for j = 1:k
    D{j} = full(A(l:l+m(j)-1, l:l+m(j)-1));
    l = l + m(j);
    if j < k
        K = full(A(l-m(j):l-1, l:l+m(j+1)-1));
        [row, col] = find(K);
        G{j} = K(min(row):end, 1:max(col));

        K = full(A(l:l+m(j+1)-1, l-m(j):l-1));
        [row, col] = find(K);
        F{j+1} = K(1:max(row), min(col):end);
    end
end
