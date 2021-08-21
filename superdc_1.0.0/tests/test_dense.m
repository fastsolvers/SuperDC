%%%%%%%% Dense matrix test %%%%%%%%

addpath(genpath('../../superdc_1.0.0'))
addpath(genpath('../../fmm1d'))
addpath(genpath('../../hss'))

fprintf('***** Dense matrix test *****\n')

%%% Toeplitz matrix example in Fourier space
fprintf(['---------------------------------------------------------\nA random Toeplitz matrix treated as DENSE\n',...
    '   with dense FFT to convert to Cauchy-like in Fourier space\n',...
    '   followed by dense HSS construction -- slow for large n\n',...
    'The compression is based on Matlab SVD with truncation\n',...
    'For nearly O(n) version, see [Xia, et al., SIMAX 33 (2012)]\n',...
    'For fast HSS construction for other matrix types\n'...
    '   and faster rank-revealing factorizations\n',...
    '   please contact the developers\n---------------------------------------------------------\n'])

n = 8192;  % matrix size
v = randn(n,1);
T = toeplitz(v); % symmetric Toeplitz matrix

fprintf('--- Dense conversion to Cauchy-like ---\n')
A = ifft((ifft(T)*sqrt(n))')*sqrt(n);

%%% The following is to make sure the matrix is "perfectly" Hermitian
if ~ishermitian(A)
    A = (A+A')/2;
    fprintf('!!! ATTN: A is not "perfectly" Hermitian and is replaced by (A+A*)/2 !!!\n');
end

%%% Dense HSS construction
m0 = 2^11;  % HSS leaf level diagonal block size
[tr,m] = npart(n,m0); % HSS partition info

tol = 1e-10; % HSS construction & deflation tolerance

fprintf('--- Dense symmetric HSS construction ---\n')
[D,U,R,B] = mat2hsssym(A,tr,m,'tol',tol);

%%% SuperDC eigensolver
accchk = 1;  % accuracy check switch

fprintf('--- SuperDC ---\n')

tic
[d,Q] = superdc(D,U,B,R,tr,tol);
tsuperdc = toc;

fprintf('SuperDC      time: %11.4e\n',tsuperdc);

%%% Application of the eigenmatrix to vectors.
x = randn(n,1);

tic
y = superdcmv(Q,x,0);
tmv = toc;

fprintf('Q mat-vec    time: %11.4e\n',tmv);

%%% SuperDC accuracy check ( decomposition residual,  loss of orthogonality )
if accchk == 1    
    % Matlab eig
    fprintf('\n--- Matlab eig ---\n')
    
    tic
    [~,D0] = eig(A);
    teig = toc;
    
    fprintf('eig          time: %11.4e\n',teig);
    
    d0 = diag(D0);
    err = norm(d-d0)/n/norm(d0);
    fprintf('\nError measurement: %11.4e\n',err)
    
    fprintf('To check residual and loss of orthogonality (slow!), comment out lines below this message...\n')
    %[res,oth] = superdcacc(A,d,Q,m0);
    %fprintf('Residual:        %11.4e\nOrthogonality:   %11.4e\n',res,oth); 
end
