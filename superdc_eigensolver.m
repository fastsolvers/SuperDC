%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  SUPERDC eigensolver interface   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lam, Q, stats] = superdc_eigensolver(type, tol, varargin)

%%% Input:
%%% tol: tolerance,

%%% type: 

%%%                   'hss':     varargin{:} = [D, U, R, B, tr] 
%%%                                         (   symmetric HSS form of A,     recommended HSS leaf block size 2048   )


%%%            'banded':     varargin{:} = A, or
%%%                                 varargin{:} = [M,d,n] where A = spdiags(M,d,n,n)


%%%         'blktridiag':     varargin{:} = [T, S]
%%%                                  T : cell array, k diagonal blocks    (   recommended block size 512   )
%%%                                  S : cell array, k-1 subdiagonal blocks


%%%             'toeplitz':     varargin{:} = t  
%%%                                 (first column/row of A)
%%%   !!!!!!!!!! ATTENTION: this option requires the Toeplitz Solver package from  http://www.mathcs.emory.edu/~yxi26/  



%%% Output:
%%% lam: eigenvalues, column vectors
%%% Q: hss structured eigenmatrix
%%% stats: performance report (flops, time, storage, etc.)



%%  HSS form input / construction

if strcmp(type, 'hss')
    [D, U, R, B, tr] = varargin{:};
    
    
elseif strcmp(type, 'banded')
     if length(varargin) == 1
        A = varargin{:};
     else
        [M, d, n] = varargin{:};
        A = spdiags(M, d, n, n);
     end
     [T, F, G] = banded2blktri(A);
     [D, U, R, B, ~, ~, tr] = blktri2hss(T, F, G, 4);
     
     
elseif strcmp(type, 'blktridiag')
     [T, S] = varargin{:};
     G = cell(length(T),1);
     F = cell(length(T),1);
     for i=1:length(T)-1
         G{i} = S{i}.';
         F{i+1} = S{i};
     end
     [D, U, R, B, ~, ~, tr] = blktri2hss(T, F, G, 4);

     
elseif strcmp(type, 'toeplitz')
     t = varargin{:};
     t = reshape(t, length(t), 1);
     n = length(t);
     t = reshape(t, n, 1);
     r = min(2048, floor(n/8));
     [tr,m] = npart(n,r); 
     
     width = 40;
     [d1, u, v, dd] = toep2cauchysym(t);
     [D, U, R, B] = Cau2hsssym(t, d1, u, v, dd, n, tr, m, 'tol', tol, width);
end





%% SuperDC eigensolver
N = 2^10;
[lam, Q, time_divide, flops_divide, time_conquer, flops_conquer, rho] = superdc(D, U, B, R, tr, tol, N);




%% performance statistics
stats = struct;

stats.total_time = time_divide + time_conquer;
stats.total_flops = flops_divide + flops_conquer;

stats.time_divide = time_divide;
stats.time_conquer = time_conquer;

stats.flops_divide = flops_divide;
stats.flops_conquer = flops_conquer;

stats.storage = Qstorage(Q);
stats.rho = rho;

