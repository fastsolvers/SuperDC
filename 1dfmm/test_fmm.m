
addpath('../1dfmm')
fun = 3;
% nx = 5000;
% ny = 5000;
% x = 1000*abs(randn(1,nx));
% y = 1000*abs(randn(1,ny));
% x = sort(x).';
% y = sort(y).' ;

n = 10000;
s = 10000*abs(randn(n,1));
s = sort(s);
y = s(1:2:end);
x = s(2:2:end);
scaling = 1;
nx = length(x);
ny = length(y);
org = [1:nx]';
dif = x-y;

k = floor(nx/2);
J = randsample(ny-1,k);
org(J) = org(J)+1;
dif(J) = x(J) - y(J+1);
q = ones(size(y,1),1);


T = 1;
switch fun
case 1
    if T == 0
        Ker = 1 ./ (x  - y.');
    else
        Ker = 1 ./ (y  - x.');
    end
case 2
    if T == 0
        Ker = 1 ./ (x  - y.').^2;
    else
        Ker = 1 ./ (y  - x.').^2;
    end
case 3
    if T == 0
        Ker = log(abs((x  - y.')));
    else
        Ker = log(abs((y  - x.')));
    end
end  
Ker(isinf(Ker)) = 0;
phi = Ker * q;



for r = 20:10:60
    fprintf('r = %d\n', r )    
    tic
%     [z] = fmm1d( r, x, y, q, fun, scaling ); 
    [zl, zu, nflops] = trifmm1dlu_local_shift_2( r, y, x, q, dif, org, fun, scaling);
    time3 = toc;
    error = norm( phi - zl - zu ) / (norm( phi ));
  
    fprintf('fmm error: |phi-z|/|phi| = %e\n',error );
    fprintf('fmm time: %e\n', time3);
    fprintf('\n\n')
end