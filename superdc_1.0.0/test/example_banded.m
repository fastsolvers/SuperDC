addpath(genpath('../../../SuperDC'))
addpath(genpath('../../../SuperDC-main'))

%%%%%%%%%%%%%%%%%%%%%% example 3 ( banded )   %%%%%%%%%%%%%%%%%%%%%%%%



%% Example set up
tol = 1e-10;
n = 2^13;       % matrix size
bw = 5;           % bandwidth
M = randn(n, 2*bw+1);
d = -bw:bw;
A = spdiags(M, d, n, n);
A = (A + A.') / 2;    





%% SuperDC eigensolver
[lam, Q, stats] = superdc_eigensolver('banded', tol, A);
t_superdc = stats.total_time;




%% Application of the eigenmatrix to vectors.
X = randn(n,10);
[Y, stats2] = superdc_matvec(Q, X, 0, 'banded');




%% MATLAB eig
Anorm = normest(A);
tic
[~, D2] = eig(full(A));
t_eig = toc;
lam2 = diag(D2);



%% SuperDC accuracy check ( decomposition residual,  loss of orthogonality )

 % 0 < sample_ratio <= 1    
sample_ratio = 0.05;  

% sample without replacement (n * sample_ratio) eigenvectors to compute residual and orthogonality
[decomposition_residual, loss_orthogonality] = superdc_accuracy(Q, lam, sample_ratio, 'banded', A);




%% Results 
fprintf('\n\n================== example 3 : banded =====================\n\n'  )
fprintf('================== matrix size: %i =====================\n\n', n)
fprintf('eigenvalues error: %e\n', norm(lam - lam2) / (n * norm(lam2)));
fprintf('decomposition residual: %e\n', max(decomposition_residual) / Anorm);
fprintf('loss of orthogonality: %e\n', max(loss_orthogonality));
fprintf('superdc time : %.1f\n', t_superdc);
fprintf('eig time : %.1f\n\n\n', t_eig);



