%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% compute vhat by Lowner's formula %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function   [v, nflops] = vhat(d, lam, tau, org, w, N)
%%% Input:
%%% d lam w: as in secular equation
%%% K: indices of dense eigenvectors
%%% Kc: indices of structured eigenvectors
%%% union(K, Kc) = 1:n
%%% QK: dense eigenvectors
%%% N: size threshold to use fmm

%%% Output
%%% v: output of Lowner's formula

if nargin < 6
    N = 2^10;
end

n = length(lam);
lam = reshape(lam, [n,1]);
d = reshape(d, [n,1]);
r = 50;
nflops = 0;

if n < N
    Dlam = d - d(org).';
    Dlam = bsxfun(@minus, Dlam, tau.');
    D = d - d.';
    D(1:(n+1):end) = ones(n,1);
    v = (log(abs(Dlam))-log(abs(D))) * ones(n,1); v1 = v;
    v = exp(1/2 * v);
    v(w < 0) = - v(w < 0);
    
    nflops = nflops + flops('prod', D, 'n', v, 'n') + numel(D) + 2*numel(v);
else            
    e1 = ones(n, 1);
    [v0, nflops1] = fmm1d_local_shift_2(r, d, lam, e1, tau, org, 3);
    nflops = nflops + nflops1;

    [vd, nflops1] = fmm1d_local_shift(r, d, d, e1, zeros(n,1), [1:n].', 3);
    nflops = nflops + nflops1;
    
    v = exp(1/2 * (v0 - vd));
    v(w < 0) = - v(w < 0);
    nflops = nflops + 4*numel(v);
end




