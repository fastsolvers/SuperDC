%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% compute columns norms of Cauchy matrix ( v_i / (d_i - lam_j) ) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [s, nflops] = colnorms(d, lam, tau, org, v, N)
%%% Input:
%%% d, lam, tau, org, v : from secular.m and rootfinder.m
%%% N: size threshold to use fmm

%%% Output
%%% s: norm of each columns of ( v_i / (d_i - lam_j) )

if nargin < 6
    N = 2^10;
end

n = length(lam);
v = reshape(v, [n, 1]);
d = reshape(d, [n, 1]);
lam = reshape(lam, [n, 1]);
r = 50;
nflops = 0;

if n < N 
    S = d(org) - d.';
    S = bsxfun(@plus, S, tau);
    S = 1 ./ S.^2;
    s = S * (v.^2);
    nflops = nflops + flops('prod', S, 'n', v, 'n');
    
    s = 1 ./ sqrt(s);
    nflops = nflops + numel(s);
    
else                     
    [s, nflops1] = fmm1d_local_shift(r, lam, d, v.^2, tau, org, 2);
    nflops = nflops + nflops1;

    s = 1 ./ sqrt(s); 
    nflops = nflops + numel(s);
end






