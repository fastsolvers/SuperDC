function [res,oth] = superdcacc(A,d,Q,m0)

% Residual

n = length(d);
na = norm(A);

for i = 1:n
    e = sparse(n,1); e(i) = 1;
    q = superdcmv(Q,e,0,m0);
    rv(i) = norm(A*q-d(i)*q)/n/na;
    qq = superdcmv(Q,q,1,m0);
    ov(i) = norm(qq-e)/n;
end

res = max(rv);
oth = max(ov);