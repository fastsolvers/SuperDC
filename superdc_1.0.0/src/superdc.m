function [lam,Q,time_divide,flops_divide,time_conquer,flops_conquer,rho] = superdc(D,U,B,R,tr,tol,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SuperDC eigensolver %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input:
% D,U,B,R: HSS generators
% tr:      HSS tree
% tol:     tolerance for deflation
% N:       (optional) size threshold to use fmm; default 1024
%
%%% Output:
% Lam:     eigenvalues
% Q:       hss structured eigenmatrix
% Other outputs for results reporting purpose

if nargin < 7
    N = 2^10;
end


%% HSS tree
m = partitionsz(D, tr);
rg = indrange(tr,m);
ch = child(tr);
k = length(tr);


%% dividing stage
t = clock;
[D, ~, Z, desc, flops_divide] = divide2(D, U, B, R, tr);
time_divide = etime(clock, t);


%% conquering stage
flops_conquer = 0;
Q0 = cell(k,1);         
Lam = cell(k,1);
rho = cell(k,1);
t = clock;
for i = 1:k
    if isempty(ch{i})
        [Q0{i}, lam] = eig(D{i});
        flops_conquer = flops_conquer + 4/3 * size(D{i}, 1)^3;
        Lam{i} = diag(lam);
    else
        c1 = ch{i}(1);
        c2 = ch{i}(2);
        
        [Z{i}, nflops1] =  superdcmv_desc(Q0, Z{i}, tr, i, 1, rg, desc, N);
        flops_conquer = flops_conquer + nflops1;
        
        Lam{i} = [Lam{c1}; Lam{c2}];
        Lam{c1} = [];
        Lam{c2} = [];
        r = size(Z{i}, 2);
        rho{i} = zeros(r, 1);
        for j = 1:r    
            [Lam{i}, Q0{i}{j}, nflops1, rho1] = secular(Lam{i}, Z{i}(:, j), tol, N);
            flops_conquer = flops_conquer + nflops1;
            rho{i}(j) = rho1;
            if j < r
                [Z{i}(:, j+1:r), nflops1] = superdcmv_cauchy(Q0{i}{j}, Z{i}(:, j+1:r), 1, N);  
                flops_conquer = flops_conquer + nflops1;
            end
        end
    end
end
time_conquer = etime(clock, t);


%% eigenvalues and eigenmatrix
lam = Lam{k};
[lam, I] = sort(lam, 'ascend');
Q = {Q0, I, tr, m};
