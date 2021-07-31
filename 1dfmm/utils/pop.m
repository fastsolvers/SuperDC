function [S,U] = pop(S)

n = length(S);
% U = full(sparse(S{n}(:,1),S{n}(:,2),S{n}(:,3)));
% U = accumarray(S{n}(:,1:2),S{n}(:,3));

U = reshape(S{n}(3:end),S{n}(1:2)');
S(n) = []; % new: 2008-9-15
% is = is(1:n);
% i = it(end);
% it = it(1:end-1);