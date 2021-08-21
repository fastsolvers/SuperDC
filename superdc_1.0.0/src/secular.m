function [Lam, Q, nflops, percent] = secular(d, v, tol, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% roots of secular equation w/ deflation %%%%%
%%%%% 1 + \sum_j (v_j^2 / (lam1_j - lam)) = 0%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Input:
% d,v: as in secular equation
% tol: tolerance for deflation
% N:   size threshold to use fmm
%
%%% Ouput:
% Lam: roots of secular equation
% Q:   structured eigen matrix

if nargin < 4
    N = 2^10;
end

nflops = 0;
n = length(v);
d = reshape(d, [n, 1]);
v = reshape(v, [n, 1]);

%%%% step 1: deflate small vi
T = find(abs(v) < tol)';
Tc = setdiff(1:n, T);
Lam1 = d(T);
Lam2 = [];
Lam3 = [];
d2 = d(Tc);
v2 = v(Tc);
n1 = length(T);
T = [T, Tc];

%%%% step 2: deflate close eigenvalue
if n1 < n
    v2c = conj(v2) ./ abs(v2);
    v2 = abs(v2);
    [d2, I] = sort(d2, 'ascend');
    v2 = v2(I);
    
    p = 1; J = []; G = [];
    for j = 2:length(d2)
        s = v2(p);
        c = v2(j);
        tau = sqrt(c^2+s^2);
        t = d2(j) - d2(p);
        c = c / tau;
        s = s / tau;
        
        if abs(t*c*s) < tol                 % close eigenvalues, deflates d2(p)
            v2(p) = 0;
            v2(j) = tau;
            t = c^2 * d2(p) + s^2 * d2(j);
            d2(j) = s^2 * d2(p) + c^2 * d2(j);
            d2(p) = t;
            G = [G; p, j, c, s];                % Givens rotations
            J = [J, p];
            p = j;
        else
            p = j;
        end
    end
    Jc = setdiff(1:n-n1, J);
    Lam2 = d2(J);
    d3 = d2(Jc);
    v3 = v2(Jc);
    
    n2 = length(J);
    n3 = n - n1 - n2;
    J = [J, Jc];
    
    if n3
        [Lam3, tau, org, nflops1, percent] = rootfinder(d3, v3, N);
        nflops = nflops + nflops1;
        
        [v3_hat, nflops1] = vhat(d3, Lam3, tau, org, v3, N);
        nflops = nflops + nflops1;
        
        [s3, nflops1] = colnorms(d3, Lam3, tau, org, v3_hat, N);
        nflops = nflops + nflops1;
    end
    
end


Lam = [Lam1; Lam2; Lam3];
if n1 < n
    Q3 = {v3_hat, s3, d3, Lam3, tau, org};
    Q = {Q3, J, G, I, v2c, T, [n, n1, n2, n3]};
else
    Q3 = {}; J = []; G = []; I = []; v2c = [];
    Q = {Q3, J, G, I, v2c, T, [n, n1, 0, 0]};
end












%%%%%%%%% test and check %%%%%%%%%

%%%%%%%% for error test %%%%%%%%%%
%  if n3
%      S = d3 - d3(org).';
%      S = bsxfun(@minus, S, tau.');
%      S = 1 ./ S;
%      S = bsxfun(@times, S, v3_hat);
%      S = bsxfun(@times, S, s3.');
%      res = max(sqrt(sum((S' * S - eye(n3)).^2, 1))) / n3;
%      fprintf('rootfinder: loss of orthogonality: %e, %i\n', res, n3);
%  end

%%%%%%%% for error test %%%%%%%%%%


% fprintf('node: %i, stage: %i, size: %i\n', [node, stage, n3]);
% A1 = diag(d3) + v3 * v3.';
% A2 = diag(d3) + v3_hat * v3_hat.';
%
% fprintf('secualrEq: max(abs((v - v_hat))) : %e \n', max(abs(v3 - v3_hat)));
% fprintf('secualrEq: min(abs((v - v_hat))) : %e \n', min(abs(v3 - v3_hat)));
% fprintf('secualrEq: max(abs((v - v_hat))) / norm(v): %e \n\n', norm(v3 - v3_hat, Inf) / (n3 * norm(v3)));
% %
% fprintf('secularEq: norm(A - A_hat) / norm(A): %e \n', norm(A1 - A2) / norm(A1));
% fprintf('secularEq: norm(A): %e\n', norm(A1));
% fprintf('secularEq: norm(A - A_hat): %e \n', norm(A1 - A2));
%
% fprintf('secularEq: norm(A - Q * Lam * Q^T): %e\n', norm(A1 - Q3 * diag(Lam3) * Q3.'))
% fprintf('secularEq: norm(A_hat - Q * Lam * Q^T): %e\n\n', norm(A2 - Q3 * diag(Lam3) * Q3.'))
% fprintf('secularEq: residual of eigen decomposition1: %e\n', max(sqrt(sum((A1 * Q3 - Q3 * diag(Lam3)).^2, 1))));
% fprintf('secularEq: residual of eigen decomposition2: %e\n\n', max(sqrt(sum((A2 * Q3 - Q3 * diag(Lam3)).^2, 1))));
%
% fprintf('secularEq: relative residual of eigen decomposition1: %e\n', max(sqrt(sum((A1 * Q3 - Q3 * diag(Lam3)).^2, 1))) / (n3 * norm(A1)));
% fprintf('secularEq: relative residual of eigen decomposition2: %e\n\n', max(sqrt(sum((A2 * Q3 - Q3 * diag(Lam3)).^2, 1))) / (n3 * norm(A2)));
%
% fprintf('secularEq: norm(Q * Q^T - I): %e\n', norm(Q3 * Q3.' - eye(n3)))
% fprintf('secularEq: node: %i, stage: %i, size: %i\n', [node, stage, n]);
% fprintf('secularEq: number of eigenvalues deflated at node %i, stage %i: [[%i, %i]] (out of %i) \n', [node, stage, n1, n2, n]);
% fprintf('secularEq: loss of orthogonality: %e\n\n\n\n\n\n', max(sqrt(sum((Q3' * Q3 - eye(n3)).^2, 1))) / n3);
