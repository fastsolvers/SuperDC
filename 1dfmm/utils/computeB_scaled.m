function [ B ] = computeB_scaled( r, a, b, dx, dy, fun, scaling)
% center: a, diam:  dx
% fun = 1, 1/(x-y)
% fun = 2, 1/(x-y)^2
% fun = 3, ln|x-y|
% fun = 4, ln(x-y), only for dim=2

ba = b-a;
S = diag( (-1).^(r-1:-1:0) );
S = S( 1:r, r:-1:1 );
DD = diag( (-1).^(0:r-1) );
B = zeros(r);

if fun == 1
    switch scaling
        case 0
            % no scaling in B
            cc = -ones(r,1)/ba;
            for k = 1:r-1
                cc(k+1) = k*cc(k)/ba;
            end
            B = zeros(r);
            for i = 1:r
                B(i,1:r-i+1) = cc(i:end);
            end
            
        case 1
            sr =  (2*pi*r)^(1/2/r)/exp(1);    
            B(1,1) = -1/ba;
            B(1,2) = -1/ba^2/(sr*2/dy);
            B(2,1) = -1/ba^2/(sr*2/dx);
            B(2,2) = 2/ba*B(2,1)/(sr*2/dy);
            % complete first row of B
            for i = 3:r
                B( 1,i ) = (i-1)/ba*B( 1,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete second row of B
            for i = 3:r-1
                B( 2,i ) = i/ba*B( 2,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete first column of B
            for i = 3:r
                B( i,1 ) = (i-1)/ba*B( i-1,1 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete second column of B
            for i = 3:r-1
                B( i,2 ) = i/ba*B( i-1,2 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete the rest
            for k = 3:r-2
                for i = 3:r-k+1
                    B( k,i ) = (k+i-2)/ba*B( k,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
                end
            end
    end
    
    B = B*DD;
    
elseif fun == 2
    switch scaling
        case 0
            eta = 1;
            etaba = 1/( eta*ba );
            cc = 1/ba^2 * ones( r,1 );
            for k = 1:r-1
                cc(k+1) = ( k+1 )*cc(k)*etaba;
            end
            C = triu( toeplitz( cc(end:-1:1) ) );
            B = C * S;
            
        case 1
            sr =  (2*pi*r)^(1/2/r)/exp(1);    
            B(1,1) = 1/ba^2;
            B(1,2) = 2/ba^3/(sr*2/dy);
            B(2,1) = 2/ba^3/(sr*2/dx);
            B(2,2) = 3/ba*B(2,1)/(sr*2/dy); 
            % complete first row of B
            for i = 3:r
                B( 1,i ) = i/ba*B( 1,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete second row of B
            for i = 3:r-1
                B( 2,i ) = (i+1)/ba*B( 2,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete first column of B
            for i = 3:r
                B( i,1 ) = i/ba*B( i-1,1 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete second column of B
            for i = 3:r-1
                B( i,2 ) = (i+1)/ba*B( i-1,2 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete the rest
            for k = 3:r-2
                for i = 3:r-k+1
                    B( k,i ) = (k+i-1)/ba*B( k,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
                end
            end
            B = B*DD;
    end
             
%     if scaling
%         B = scaleb(B,r,dx,dy);
%     end
elseif fun == 3
    switch scaling
        case 0
            eta = 1;
            etaba = 1/( eta*ba );
            cc = -etaba * ones( r,1 ); 
            for k = 2:r-1
                cc(k+1) = (k-1)*cc(k)*etaba;
            end
            cc(1) = log( abs(ba) );
            C = triu( toeplitz( cc(end:-1:1) ) );
            B = C * S;
        case 1
           sr =  (2*pi*r)^(1/2/r)/exp(1);
            B(1,1) = log(abs(a-b));
            B(1,2) = -1/ba/(sr*2/dy);
            B(2,1) = -1/ba/(sr*2/dx);
            B(2,2) = 1/ba*B(2,1)/(sr*2/dy);
            % complete first row of B
            for i = 3:r
                B( 1,i ) = (i-2)/ba*B( 1,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete second row of B
            for i = 3:r-1
                B( 2,i ) = (i-1)/ba*B( 2,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete first column of B
            for i = 3:r
                B( i,1 ) = (i-2)/ba*B( i-1,1 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete second column of B
            for i = 3:r-1
                B( i,2 ) = (i-1)/ba*B( i-1,2 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete the rest
            for k = 3:r-2
                for i = 3:r-k+1
                    B( k,i ) = (k+i-3)/ba*B( k,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
                end
            end
            B = B*DD;
    end
    
%     if scaling
%         B = scaleb(B,r,dx,dy);
%     end

elseif fun == 4
    switch scaling
        case 0
            eta = 1;
            etaba = 1/( eta*ba );
            cc = -etaba * ones( r,1 ); 
            for k = 2:r-1
                cc(k+1) = (k-1)*cc(k)*etaba;
            end
            cc(1) = log(ba);
            C = triu( toeplitz( cc(end:-1:1) ) );
            B = C * S;
        case 1
           sr =  (2*pi*r)^(1/2/r)/exp(1);
            B(1,1) = log(a-b);
            B(1,2) = -1/ba/(sr*2/dy);
            B(2,1) = -1/ba/(sr*2/dx);
            B(2,2) = 1/ba*B(2,1)/(sr*2/dy);
            % complete first row of B
            for i = 3:r
                B( 1,i ) = (i-2)/ba*B( 1,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete second row of B
            for i = 3:r-1
                B( 2,i ) = (i-1)/ba*B( 2,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete first column of B
            for i = 3:r
                B( i,1 ) = (i-2)/ba*B( i-1,1 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete second column of B
            for i = 3:r-1
                B( i,2 ) = (i-1)/ba*B( i-1,2 )/(sr*2/dx)/(i-1)*(1-1/(i-1))^(i-2);
            end
            % complete the rest
            for k = 3:r-2
                for i = 3:r-k+1
                    B( k,i ) = (k+i-3)/ba*B( k,i-1 )/(sr*2/dy)/(i-1)*(1-1/(i-1))^(i-2);
                end
            end
            B = B*DD;
    end

end
