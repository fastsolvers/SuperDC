function nflops = flops(type,A,transA,B,transB,varargin)

if strcmp(type,'prod')
    if lower(transA) == 'n'
        m = size(A,1); n = size(A,2);
    else
        m = size(A,2); n = size(A,1);
    end
    if lower(transB) == 'n'
        p = size(B,2);
    else
        p = size(B,1);
    end
    nflops = double(m*(2*n-1)*p);
    if nargin > 5 % multiple/nested products
        for i = 1:length(varargin)/2
            m = p;
            if lower(varargin{2*i}) == 'n'
                n = size(varargin{2*i-1},1); p = size(varargin{2*i-1},2);
            else
                n = size(varargin{2*i-1},2); p = size(varargin{2*i-1},1);
            end
            nflops = nflops+double(m*(2*n-1)*p);
        end
    end
    return
end
if strcmp(type,'prodsym')
    [m,n] = size(A);
    if lower(transA) == 'n' % A*A'
        nflops = double((2*n-1)*m*(m+1)/2);
    else
        nflops = double((2*m-1)*n*(n+1)/2);
    end
    return
end
if strcmp(type,'sum')
    nflops = double(numel(A));
    return
end
if strcmp(type,'sumsym')
    nflops = double(numel(A)/2+size(A,1)/2);
    return
end
if strcmp(type,'mv')
    % Ax, note; x is never included as an argument
    if nargin > 3 % multiple/nested products, multiplied from right to left
        nflops = 0;
        varargin = [{A} {transA} {B} {transB} varargin];
        for i = length(varargin)/2:-1:1 % the last argument (x) is skipped
            if lower(varargin{2*i}) == 'n'
                m = size(varargin{2*i-1},1); n = size(varargin{2*i-1},2);
            else
                m = size(varargin{2*i-1},2); n = size(varargin{2*i-1},1);
            end
            nflops = nflops+double(m*(2*n-1));
        end
    else
        if lower(transA) == 'n'
            m = size(A,1); n = size(A,2);
        else
            m = size(A,2); n = size(A,1);
        end
        nflops = double(m*(2*n-1));
    end
    return
end
if strcmp(type,'rdiv')
    if lower(transA) == 'n'
        if strcmp(transB,'tri')
            nflops = double(size(A,1)*size(B,1)^2);
        else
            nflops = double(size(A,1)*2/3*size(B,1)^3);
        end
    else
        if strcmp(transB,'tri')
            nflops = double(size(A,2)*size(B,1)^2);
        else
            nflops = double(2*size(A,2)*size(B,1)^2+2/3*size(B,1)^3);
        end
    end
    return
end

if strcmp(type,'ldiv')
    if lower(transB) == 'n'
        if strcmp(transA,'tri')
            nflops = double(size(B,2)*size(A,1)^2);
        else
            nflops = double(size(B,2)*2/3*size(A,1)^3);
        end
    else
        if strcmp(transA,'tri')
            nflops = double(size(B,1)*size(A,1)^2);
        else
            nflops = double(2*size(B,1)*size(A,1)^2+2/3*size(A,1)^3);
        end
    end
    return
end

if strcmp(type,'chol')
    nflops = double(size(A,1)^3/3)-size(A,1)^2+size(A,1)*2/3;
    return
end

if strcmp(type,'lu')
    nflops = 2*double(size(A,1)^3/3)-size(A,1)^2/2-size(A,1)/6;
    return
end


if strcmp(type,'qr')
    nflops = 2*size(A,1)*size(A,2)^2 - 2/3*size(A,2)^3;
    return
end

if strcmp(type,'ldl')
    nflops = size(A,1)^3/3-size(A,1)^2+size(A,1)*2/3;
    return
end
