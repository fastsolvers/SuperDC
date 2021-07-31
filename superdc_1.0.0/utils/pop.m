function [S,U] = pop(S)
%%% pop out the top element of stack S
n = length(S);
U = reshape(S{n}(3:end),S{n}(1:2)');
S(n) = []; 