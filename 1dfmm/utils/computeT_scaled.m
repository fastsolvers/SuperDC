function [ T ] = computeT_scaled( r, a, b, dx, dy, scaling)
% compute translation matrix: a->b
% center of cluster x: a
% diam of cluster x:  dx
% init and compute diagonal of T, simultaneously

sr =  (2*pi*r)^(0.5/r)/exp(1);

T = diag( (dx/dy).^(0:r-1) );
ab = a-b;
T(1,2) = ab*(sr*2/dy);

k = 1;i = 2;
if abs(T(k,i)-(i-1)^(i-1)/(k-1)^(k-1)*sr^(i-k)*(dx/2)^(k-1)/(dy/2)^(i-1)*ab^(i-k)/factorial(i-k)) > 0
    (i-1)^(i-1)/(k-1)^(k-1)
    dx^(k-1)
    error
end

%scaling = 0;

switch scaling
    case 1
        % complete first row of T
        for i = 3:r
            T( 1,i ) = ab/(i-1)*T( 1,i-1 )*(sr*2/dy)*(i-1)/(1-1/(i-1))^(i-2);
            k = 1;
            %T(k,i) = (i-1)^(i-1)/(k-1)^(k-1)*sr^(i-k)*dx^(k-1)/dy^(i-1)*ab^(i-k)/factorial(i-k);
            if abs(T(k,i)) > 1
               T(k,i)
                abs(T(k,i)-(i-1)^(i-1)/(k-1)^(k-1)*sr^(i-k)*(dx/2)^(k-1)/(dy/2)^(i-1)*ab^(i-k)/factorial(i-k))
                %T(k,1:i)
                %[T(k,i) (i-1)^(i-1)/(k-1)^(k-1)*sr^(i-k)*(dx/2)^(k-1)/(dy/2)^(i-1)*ab^(i-k)/factorial(i-k)]
                error
            end

            % L = abs(ab)/(i-1)*abs(T( 1,i-1 ))*(sr*2/dy)*(i-1)/(1-1/(i-1))^(i-2);
% theta = (i-1)*angle(ab);
% T( 1,i ) = L*complex(cos(theta),sin(theta));
        end
        % complete the rest
        for k = 2:r
            for i = k+1:r
                T( k,i ) = ab/(i-k)*T( k,i-1 )*(sr*2/dy)*(i-1)/(1-1/(i-1))^(i-2);
%     L = abs(ab)/(i-k)*abs(T( k,i-1 ))*(sr*2/dy)*(i-1)/(1-1/(i-1))^(i-2);
%     theta = (i-k)*angle(ab);
%     T( k,i ) = L*complex(cos(theta),sin(theta));    
%             if abs(T(k,i)) > 1
%                T(k,i),k,i,sr,dx,dy,ab
%                 abs(T(k,i)-(i-1)^(i-1)/(k-1)^(k-1)*sr^(i-k)*(dx/2)^(k-1)/(dy/2)^(i-1)*ab^(i-k)/factorial(i-k))
%                 %T(k,1:i)
%                 %[T(k,i) (i-1)^(i-1)/(k-1)^(k-1)*sr^(i-k)*(dx/2)^(k-1)/(dy/2)^(i-1)*ab^(i-k)/factorial(i-k)]
%                 error
%             end
            %T(k,i) = (i-1)^(i-1)/(k-1)^(k-1)*sr^(i-k)*dx^(k-1)/dy^(i-1)*ab^(i-k)/factorial(i-k);
            end
        end


    case 0
        % no scaling in T
        cc = ones(r,1);
        for k = 1:r-1
            cc(k+1) = ab*cc(k)/k;
        end
        T = zeros(r);
        for i = 1:r
            T(i,i:r) = cc(1:r-i+1);
        end
        
    case 2
        % straightfoward calculation to show bound of T        
        ma = sr/dx*2;
        Sa = diag([ 1, (ma*(1:r-1)).^(1:r-1) ]);
        mb = sr/dy*2;
        Sb = diag([ 1, (mb*(1:r-1)).^(1:r-1) ]);
        
        cc = ones(r,1);
        for k = 1:r-1
            cc(k+1) = ab*cc(k)/k;
        end
        T = zeros(r);
        for i = 1:r
            T(i,i:r) = cc(1:r-i+1);
        end
        
        T = 1./Sa*T*Sb;
        
        
end
