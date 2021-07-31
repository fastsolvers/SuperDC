addpath(genpath('../../../SuperDC'))
addpath(genpath('../../../SuperDC-main'))

%%%%%%%%%%%%%%%%%%%%%% example 5 ( toeplitz )   %%%%%%%%%%%%%%%%%%%%%%%%


%% Example set up
tol = 1e-10;
n = 2^11;

% % Prolate matrix
% a = 1/4;
% t = [2*a, sin(2*a*[1:n-1]*pi) ./ ([1:n-1]*pi)].';
% A = toeplitz(t);
% load prolate_hss.mat
% [D, U, R, B, tr] = prolate_hss{:};      % HSS form of the Cauchy-like matrix after transformation, recommended HSS leaf block size = 2048


% % KMS matrix
rho = 1/2;
t = rho.^([0:n-1].');
A = toeplitz(t);
load kms_hss.mat
[D, U, R, B, tr] = kms_hss{:};      % HSS form of the Cauchy-like matrix after transformation, recommended HSS leaf block size = 2048





%% SuperDC eigensolver
[lam, Q, stats] = superdc_eigensolver('hss', tol, D, U, R, B, tr);
t_superdc = stats.total_time;

% %  !!!!!!!!!! ATTENTION: this option requires the Toeplitz Solver package from  http://www.mathcs.emory.edu/~yxi26/  

% [lam, Q, stats] = superdc_eigensolver('toeplitz', tol, t); 
% t_superdc = stats.total_time;






%% Application of the eigenmatrix to vectors.
X = randn(n,10);
[Y, stats2] = superdc_matvec(Q, X, 0, 'toeplitz');




%% MATLAB eig
Anorm = normest(A);
tic
[~, D2] = eig(A);
t_eig = toc;
lam2 = diag(D2);




%% SuperDC accuracy check ( decomposition residual,  loss of orthogonality )

 % 0 < sample_ratio <= 1
sample_ratio = 0.05; 

% sample without replacement (n * sample_ratio) eigenvectors to compute residual and orthogonality
[decomposition_residual, loss_orthogonality] = superdc_accuracy(Q, lam, sample_ratio, 'toeplitz', t);





%% Results 

fprintf('\n\n================== example 5 : toeplitz =====================\n\n'  )
fprintf('================== matrix size: %i =====================\n\n', n)
fprintf('eigenvalues error: %e\n', norm(lam - lam2) / (n * norm(lam2)));
fprintf('decomposition residual: %e\n', max(decomposition_residual) / Anorm);
fprintf('loss of orthogonality: %e\n', max(loss_orthogonality));
fprintf('superdc time : %.1f\n', t_superdc);
fprintf('eig time : %.1f\n\n\n', t_eig);

