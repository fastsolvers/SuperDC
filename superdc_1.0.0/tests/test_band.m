%%%%%%%%%%%%%%%%%%%%%% Banded matrix test %%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../../superdc_1.0.0'))
addpath(genpath('../../fmm1d'))

%%% Banded matrix example
n = 4096;  % matrix size
w = 1;     % half bandwidth
A = 3*eye(n)-diag(ones(n-1,1),1)-diag(ones(n-1,1),-1);

m0 = 2^11;  % HSS leaf level diagonal block size
[tr,m] = npart(n,m0); % HSS partition info

%%% band2hss
[D,U,R,B] = band2hss(A,w,tr,m);

%%% SuperDC eigensolver
tol = 1e-10;
accchk = 1;

fprintf('SuperDC ...\n')
tic
[d,Q,time_divide,flops_divide,time_conquer,flops_conquer,rho] = superdc(D,U,B,R,tr,tol);
tsuperdc = toc;
fprintf('SuperDC     time: %11.4e\n',tsuperdc);

%%% Application of the eigenmatrix to vectors.
x = randn(n,1);
tic
[y,nflops] = superdcmv(Q,x,0);
tmv = toc;
fprintf('Q mat-vec   time: %11.4e\n',tmv);

%%% SuperDC accuracy check ( decomposition residual,  loss of orthogonality )
if accchk == 1    
    % MATLAB eig
    fprintf('\nMatlab eig ...\n')
    tic
    [~,D0] = eig(A);
    teig = toc;
    fprintf('eig         time: %11.4e\n',teig);
    
    d0 = diag(D0);
    err = norm(d-d0)/n/norm(d0);
    fprintf('\nEigenvalue error: %11.4e\n',err)
    
    fprintf('To check residual and loss of orthogonality (slow!), comment out lines below this message...\n')
    %[res,oth] = superdcacc(A,d,Q,m0);
    %fprintf('Residual:        %11.4e\nOrthogonality:   %11.4e\n',res,oth); 
end