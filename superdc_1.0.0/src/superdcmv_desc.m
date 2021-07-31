%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   mat-vec of eigenmatrix of descendants of i  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X, nflops] = superdcmv_desc(Q, X, tr, i, t, rg, desc, N)
%%% Input:
%%% Q: hss structured cauchylike eigenmatrix of descendants of i
%%% X: vectors to be multiplied
%%% t = 0: not transpose
%%% t = 1: transpose
%%% N: size threshold to use fmm
%%% tr: hss tree
%%% rg: indices range of each hss leaf block 
%%% desc: smallest descendant of node i in the hss tree tr

%%% Output
%%% t = 0: X = Q * X 
%%% t = 1: X = Q^T * X;

if nargin < 8
    N = 2^10;
end

rg = rg - rg(i,1) + 1;
k = desc(i);           
nflops = 0;

if t == 0        
    for j = i-1:-1:k
        [X(rg(j, 1):rg(j, 2), :), nflops1] = superdcmv_node(Q{j}, X(rg(j, 1):rg(j, 2), :), tr, j, 0, N);
        nflops = nflops + nflops1;
    end
    
else              
    for j = k:i-1
        [X(rg(j, 1):rg(j, 2), :), nflops1] = superdcmv_node(Q{j}, X(rg(j, 1):rg(j, 2), :), tr, j, 1, N);
        nflops = nflops + nflops1;
    end
end
