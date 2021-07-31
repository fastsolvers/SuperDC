function [ U ] = computeU_scaled( x, r, a, dx, scaling)
% a: center of cluster x
% dx: diam of cluster x
% x is a column vector

if size( x,1 ) < size( x,2 )
    x = x.';
end
n = length(x);
U = ones( n, r );
sr =  (2*pi*r)^(1/2/r)/exp(1);

% compute second column of U
U(:,2) = sr*2/dx*(x-a);

%scaling = 0;

switch scaling
    case 1
        % compute the rest columns
        for k = 2:r-1
            U(:,k+1) = (1+1/(k-1))^(k-1)*sr*2/dx*(x-a).* U(:,k);
        end
    case 0
        % no scaling in U
        U = ones( n, r );
        for k = 1:r-1
            U(:,k+1) = 1/k*(x-a).* U(:,k);
        end
end

