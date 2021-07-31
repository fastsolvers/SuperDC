function S = push(S,U)

n = length(S);
% [ii,jj,v] = find(U);
% ii = [ii; size(U,1)]; jj = [jj; size(U,2)]; v = [v; 0];

% if n == 1
%     is = [0; length(ii)];
% else
%     is = [is; is(end)+length(ii)];
% end
%S{n} = [ii jj v];
S{n+1} = [size(U,1); size(U,2); reshape(U,[numel(U),1])];
% it = [it; r];

% if n == 1
%     is = [0; m];
%     S(1:m,1:m) = U;
% else
%     is = [is; is(end)+m];
% %     S = mexassign(S,U,is(n)+1,1);
%     S(is(n)+1:is(end),1:m) = U;
% end
% it = [it; r];
