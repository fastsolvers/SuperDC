addpath(genpath('../../../SuperDC'))
addpath(genpath('../../../SuperDC-main'))

%%%%%%%%%%%%%%%%%%%%%% example 4 ( block tridiagonal )   %%%%%%%%%%%%%%%%%%%%%%%%



%% Example set up
tol = 1e-10;
k = 16;                 % number of diagonal blocks
b = 512;               % diagonal blocks size   ( recommended diagonal block size 512 )
r = 5;                    % subdiagonal blocks size  
n = b * k;              % matrix size



% random symmetric block tridiagonal matrix
T = cell(k,1);        
S = cell(k-1,1);     
for i = 1:k
    Tmp = randn(b);                         % symmetric diagonal blocks
    T{i} = (Tmp + Tmp.') / 2;              
    
    if i < k
        S{i} = randn(r);                          % subdiagonal blocks
    end
end




%% SuperDC eigensolver
[lam, Q, stats] = superdc_eigensolver('blktridiag', tol, T, S);
t_superdc = stats.total_time;




%% Application of the eigenmatrix to vectors.
X = randn(n,10);
[Y, stats2] = superdc_matvec(Q, X, 0, 'blktridiag');





%% MATLAB eig
G = cell(length(T),1);
F = cell(length(T),1);
for i=1:length(T)-1
     G{i} = S{i}.';
     F{i+1} = S{i};
end
[D, U, R, B, ~, ~, tr, m] = blktri2hss(T, F, G, 4);
A = hss2symmat(D, U, R, B, tr);
Anorm = normest(A);

tic
[~, D2] = eig(A);
t_eig = toc;
lam2 = diag(D2);




%% SuperDC accuracy check ( decomposition residual,  loss of orthogonality )

% 0 < sample_ratio <= 1
sample_ratio = 0.05;        

% sample without replacement (n * sample_ratio) eigenvectors to compute residual and orthogonality
[decomposition_residual, loss_orthogonality] = superdc_accuracy(Q, lam, sample_ratio, 'blktridiag', T, S);






%% Results 
fprintf('\n\n================== example 4: block tridiagonal =====================\n\n'  )
fprintf('================== matrix size: %i =====================\n\n', n)
fprintf('eigenvalues error: %e\n', norm(lam - lam2) / (n * norm(lam2)));
fprintf('decomposition residual: %e\n', max(decomposition_residual) / Anorm);
fprintf('loss of orthogonality: %e\n', max(loss_orthogonality));
fprintf('superdc time : %.1f\n', t_superdc);
fprintf('eig time : %.1f\n\n\n', t_eig);


