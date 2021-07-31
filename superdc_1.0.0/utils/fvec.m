function [X, nflops] = fvec(b, t)


N = size(b,1);

omega = exp(pi*1i/N);

d0 = omega.^([1:N]');

X = b;

if t == 0
    X = [X(end,:); X(1:end-1,:)]; 
    X = fft(X) / sqrt(N);
    X = bsxfun(@times, X, reshape(d0,N,1));
else
    X = bsxfun(@times, X, reshape(conj(d0),N,1));
    X = ifft(X) * sqrt(N);
    X = [X(2:end,:); X(1,:)];
end

nflops = 5*N*log(N)*size(b,2)+ numel(X);