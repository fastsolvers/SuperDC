addpath(genpath('../../../SuperDC'))
addpath(genpath('../../../SuperDC-main'))

%%%%%%%%%%%%%%%%%%%%%% example 2 ( hss )   %%%%%%%%%%%%%%%%%%%%%%%%



%% Example set up

% Please provide the HSS representation of A,
% recommended HSS leaf block size = 2048
% an example is as follows

tol = 1e-6;
load kernel_hss
[D, U, R, B, tr] = kernel_hss{:};




%% SuperDC eigensolver
[lam, Q, stats] = superdc_eigensolver('hss', tol, D, U, R, B, tr);
t_superdc = stats.total_time;




%% Application of the eigenmatrix to vectors.
n = size(lam,1);
X = randn(n, 10);
[Y, stats2] = superdc_matvec(Q, X, 0, 'hss');




%% MATLAB eig
A = kernelmatrix(n);
Anorm = normest(A);
tic
[~, D2] = eig(A);
t_eig = toc;
lam2 = diag(D2);




%% SuperDC accuracy check ( decomposition residual,  loss of orthogonality )

% 0 < sample_ratio <= 1
sample_ratio = 0.05;        

% sample without replacement (n * sample_ratio) eigenvectors to compute residual and orthogonality

% [decomposition_residual, loss_orthogonality] = superdc_accuracy(Q, lam, sample_ratio, 'hss', D, U, R, B, tr);
[decomposition_residual, loss_orthogonality] = superdc_accuracy(Q, lam, sample_ratio, 'full', A);



%% Results 
fprintf('\n\n================== example 2 : hss =====================\n\n'  )
fprintf('================== matrix size: %i =====================\n\n', n)
fprintf('eigenvalues error: %e\n', norm(lam - lam2) / (n * norm(lam2)));
fprintf('decomposition residual: %e\n', max(decomposition_residual) / Anorm);
fprintf('loss of orthogonality: %e\n', max(loss_orthogonality));
fprintf('superdc time : %.1f\n', t_superdc);
fprintf('eig time : %.1f\n\n\n', t_eig);



%% generate kernel matrix
function A = kernelmatrix(n)
    x = cos(pi*(2*[0:n-1] + 1) / (2*n));
    A = sqrt(abs(x.' - x));
end
