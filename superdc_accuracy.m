%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  SUPERDC accuracy check  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [decomposition_residual, loss_orthogonality] = superdc_accuracy(Q, lam, sample_ratio, type, varargin)

%%% Input:
%%% Q: structured eigenmatrix returned by superdc_eigensolver.m
%%% lam: eigenvalues, column vectors
%%% sample_ratio: sample without replacement q = n * sample_ratio eigenvectors to compute residual and orthogonality


%%% type: 

%%%                    'hss':     varargin{:} = [D, U, R, B, tr] 
%%%                                         (symmetric HSS form of A)


%%%             'toeplitz':     varargin{:} = t   
%%%                                         (first column/row of A)


%%%            'banded':     varargin{:} = A, or
%%%                                 varargin{:} = [M,d,n] where A = spdiags(M,d,n,n)


%%%         'blktridiag':     varargin{:} = [T, S]
%%%                                  T : cell array, k diagonal blocks    (   recommended block size 512   )
%%%                                  S : cell array, k-1 subdiagonal blocks


%%%            'full':     varargin{:} = A



%%% Output:
%%% decomposition_residual : || A * q - lambda * q || for sampled eigenpairs (q, lambda)'s
%%% loss of orthogonality : || Qs^T * q - e_q || for sampled eigenvector q's


if sample_ratio > 1 || sample_ratio <=0
    sample_ratio = 1;
end

n = sum(Q{4});
q = ceil(n*sample_ratio);         
S = randsample(n,q);
S = sort(S, 'ascend').';
decomposition_residual = zeros(q,1);
loss_orthogonality = zeros(q,1);


batch_sz = ceil(min(n / 10, 2^13));
batch_ct = ceil(q / batch_sz);


for j = 1:batch_ct
    batch = (j-1)*batch_sz+1 : min(q, j*batch_sz);
    r = length(batch);
    
    Z = double([1:n]' == S(batch));
    Qz = superdc_matvec(Q, Z, 0, type);
    AQz = matvec(Qz, type, varargin{:});
    QzD = bsxfun(@times,Qz, lam(S(batch))');
    
    loss_orthogonality(batch) = sqrt(sum((Qz'*Qz - eye(r)).^2, 1)) / n;
    decomposition_residual(batch) = sqrt(sum((AQz - QzD).^2, 1)) / n;
    
    clear Z Qz AQz QzD
end










% multiply A with X
function Y = matvec(X, type, varargin)

if strcmp(type, 'hss')
    [D, U, R, B, tr] = varargin{:};
    Y = symhssmatvec2(D, U, R, B, tr, X);


elseif strcmp(type, 'toeplitz')
    t = varargin{:};
    t = reshape(t, length(t), 1);
    Y = toepmv(t, t, X);


elseif strcmp(type, 'banded')
     if length(varargin) == 1
        A = varargin{:};
     else
        [M, d, n] = varargin{:};
        A = spdiags(M, d, n, n);
     end
     Y = A * X;

elseif strcmp(type, 'blktridiag')
     [T, S] = varargin{:};
     G = cell(length(T),1);
     F = cell(length(T),1);
     for i=1:length(T)-1
         G{i} = S{i}.';
         F{i+1} = S{i};
     end
     [D, U, R, B, ~, ~, tr] = blktri2hss(T, F, G, 4);
     Y = symhssmatvec2(D, U, R, B, tr, X);


elseif strcmp(type, 'full')
    A = varargin{:};
    Y = A * X;
end
    
