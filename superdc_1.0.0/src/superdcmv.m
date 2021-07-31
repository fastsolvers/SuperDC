%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SuperDC matrix vector product %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X, nflops] = superdcmv(Q, X, t, N)
%%% Input:
%%% Q: hss sturctured eigenmatrix
%%% X: vectors to be multiplied
%%% t = 0: not transpose
%%% t = 1: transpose
%%% N: size threshold to use fmm


%%% Output
%%% t = 0: X = Q * X 
%%% t = 1: X = Q^T * X;

if nargin < 4
    N = 2^10;
end

[Q0, I, tr, m] = Q{:};
k = length(tr);
rg = indrange(tr, m);
desc = treedesc(tr);
nflops = 0;

if t == 0
    X(I, :) = X;
    [X, nflops1] = superdcmv_node(Q0{k}, X, tr, k, t, N);
    [X, nflops2] = superdcmv_desc(Q0, X, tr, k, t, rg, desc, N);
    nflops = nflops + nflops1 + nflops2;
else
    [X, nflops1] = superdcmv_desc(Q0, X, tr, k, t, rg, desc, N);
    [X, nflops2] = superdcmv_node(Q0{k}, X, tr, k, t, N);
    X = X(I, :);
    nflops = nflops + nflops1 + nflops2;
end
