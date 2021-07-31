%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  structured eigen cauchylike matrix vector product  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X, nflops] = superdcmv_cauchy(Q, X, t, N)
%%% Input:
%%% Q: hss sturctured cauchylike eigenmatrix
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


Qc = Q{1};
J = Q{2};
G = Q{3};
I = Q{4};
v2c = Q{5};
T = Q{6};
n = Q{7}(1); n1= Q{7}(2); n2= Q{7}(3); n3 = Q{7}(4); 
nflops = 0;


if t == 0
    % orthogonal Cauchy eigenmatrix
    [X(n1+n2+1:n, :), nflops1] = cauchylikematvec(Qc, X(n1+n2+1:n, :), t, N);
    nflops = nflops + nflops1;
                            
    % 2nd deflation permutation 
    X(J+n1, :) = X(n1+1:end, :);
    
    % Givens rotations
    for l = size(G,1):-1:1
        p = G(l, 1);
        j = G(l, 2);
        c = G(l, 3);
        s = G(l, 4);
        X([p+n1, j+n1], :) = [c, s; -s, c] * X([p+n1, j+n1], :);
    end
    
    % eigenvalue sorting permutation
    X(I+n1, :) = X(n1+1:end, :);
    
    % conjugate normalizer
    X(n1+1:end, :) = bsxfun(@times, X(n1+1:end, :), conj(v2c));
    nflops = nflops + (n-n1)*size(X,2);
    
    % 1st deflation permutation
    X(T, :) = X;
    
    
else
    % 1st deflation permutation
    X = X(T, :);
    
    % conjugate normalizer
    X(n1+1:end, :) = bsxfun(@times, X(n1+1:end, :), v2c);
    nflops = nflops + (n-n1)*size(X,2);

    % eigenvalue sorting permutation
    X(n1+1:end, :) = X(I+n1, :);
    
    % Givens rotations
    for l = 1:size(G,1)
        p = G(l, 1); 
        j = G(l, 2);
        c = G(l, 3);
        s = G(l, 4);
        X([p+n1,j+n1], :) = [c, -s; s, c] * X([p+n1,j+n1], :);
    end
    
    % 2nd deflation permutation
    X(n1+1:end, :) = X(J+n1, :);
    
    % orthogonal Cauchy eigenmatrix
    [X(n1+n2+1:n, :), nflops1] = cauchylikematvec(Qc, X(n1+n2+1:n, :), t, N);                                 
    nflops = nflops + nflops1;
end

