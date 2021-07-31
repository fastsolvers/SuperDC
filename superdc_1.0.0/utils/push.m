function S = push(S,U)
%%% push an element on the top of stack S

n = length(S);
S{n+1} = [size(U,1); size(U,2); reshape(U,[numel(U),1])];

