function [x,nflops]=toepmv(cv,rv,b)

nflops=0;

cv(1)=rv(1);

n =length(cv);
sz=2*n-1;

sz=length(rv);

cv = [cv;rv(end:-1:2,1)];

n1=size(b,2);

b1=[b;zeros((sz-1),n1)];

F=fft(b1);
nflops= nflops+5*sz*log(sz)*size(b1,2)/log(2);

cv=fft(cv);
for i=1:size(F,1)
    F(i,:)=cv(i)*F(i,:);
end

x=ifft(F);
nflops= nflops+5*sz*log(sz)*size(b1,2)/log(2);
x=x(1:sz,:);

end

