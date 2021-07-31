%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  SUPERDC eigenmatrix vector product  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X, stats] = superdc_matvec(Q, X, t, type)

%%% Input:
%%% Q: structured eigenmatrix returned by superdc_eigensolver.m
%%% X: vectors to be multiplied 

%%% type:  
%%%        'toeplitz', 
%%%        'hss',
%%%        'banded',
%%%        'blktridiag'


%%% Output
%%% t = 0: X = Q * X 
%%% t = 1: X = Q^T * X;
%%% stats: performance report (flops, time, storage, etc.)


N = 2^10;
stats = struct;
p = size(X,2);

time = clock;
if  nargin < 4 || ~strcmp(type, 'toeplitz')
    [X, nflops] = superdcmv(Q, X, t, N);

else
    if t == 0
        [X, nflops] = superdcmv(Q, X, 0, N);
        [X, nflops1] = fvec(X, 1);         % inverse Fourier transform
        nflops = nflops + nflops1;
    else
        [X, nflops] = fvec(X, 0);           % Fourier transform
        [X, nflops1] = superdcmv(Q, X, 1, N);
        nflops = nflops + nflops1;
    end
end
time = etime(clock, time);

stats.flops_per_matvec = nflops / p;
stats.time_per_matvec = time / p;