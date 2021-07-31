%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  mat-vec of eigenmatrix of node i  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X, nflops] = superdcmv_node(Qi, X, tr, i, t, N)
%%% Input:
%%% Qi: hss sturctured cauchylike eigenmatrix of descendants of i
%%% X: vectors to be multiplied
%%% t = 0: not transpose
%%% t = 1: transpose
%%% N: size threshold to use fmm
%%% tr: hss tree

%%% Output
%%% t = 0: X = Qi * X 
%%% t = 1: X = Qi^T * X;

if nargin < 6
    N = 2^10;
end


ch = child(tr);
nflops = 0;

if isempty(ch{i}) 
    if t == 0
        X = Qi * X;
    else
        X = Qi' * X;
    end
    nflops = flops('prod', Qi, 'n', X, 'n');
    
else                         
        r = length(Qi);
        if t == 0        
            for j = r:-1:1
                [X, nflops1] = superdcmv_cauchy(Qi{j}, X, t, N);
                nflops = nflops + nflops1;
            end
        else             
            for j = 1:r
                [X, nflops1] = superdcmv_cauchy(Qi{j}, X, t, N);
                nflops = nflops + nflops1;
            end
        end
end

