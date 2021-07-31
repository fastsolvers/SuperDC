%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  compute cauchy-like mat-vec %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Y, nflops] = cauchylikematvec(Q, X, t, N)
%%% Input:
%%% Q: {v,s,d,lam} stores the cauchylike matrix  v_i * s_j / (d_i - lam_j) 
%%% X: vectors to be multiplied
%%% t = 0: not transpose
%%% t = 1: transpose
%%% N: size threshold to use fmm

%%% Output
%%% t = 0: Y = Q * X 
%%% t = 1: Y = Q^T * X;

if nargin < 4
    N = 2^10;
end
v = Q{1};
s = Q{2};
d = Q{3};
lam = Q{4};
tau = Q{5};
org = Q{6};
n = length(v);
r = 50;
nflops = 0;

if n < N       
    switch t
        case 0 
            S = d - d(org).';
            S = bsxfun(@minus, S, tau.');
            S = 1 ./ S;
            S = bsxfun(@times, S, v);
            S = bsxfun(@times, S, s.');
            Y = S * X;
            nflops = nflops + flops('prod', S, 'n', X, 'n') + 2 * numel(S);
        case 1
            S = d(org) - d.';
            S = bsxfun(@plus, S, tau);
            S = -1 ./ S;
            S = bsxfun(@times, S, s);
            S = bsxfun(@times, S, v.');
            Y = S * X;
            nflops = nflops + flops('prod', S, 'n', X, 'n') + 2 * numel(S);
    end
    

else       
    switch t
        case 0
            X = bsxfun(@times, X, s);
            [Y, nflops1] = fmm1d_local_shift_2(r, d, lam, X, tau, org, 1);
            Y = bsxfun(@times, Y, v);
            nflops = nflops + nflops1 + numel(X) + 2*numel(Y);

        case 1
            X = bsxfun(@times, X, v);
            [Y, nflops1] = fmm1d_local_shift(r, lam, d, X, tau, org, 1);
            Y = bsxfun(@times, -Y, s);
            nflops = nflops + nflops1 + numel(X) + numel(Y);
    end
   
end





    
