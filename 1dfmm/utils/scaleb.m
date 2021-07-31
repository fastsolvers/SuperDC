function B = scaleb(B,r,dx,dy)

sr =  (2*pi*r)^(1/2/r)/exp(1);
st = sr*2/dx;
for i = 2:size(B,1)
    B(i,:) = B(i,:)/st; %B(i,:)/((i-1)*sr*2/dx)^(i-1);
    st = st*sr*2/dx*(1+1/(i-1))^(i-1)*i;
end
st = sr*2/dy;
for i = 2:size(B,2)
    B(:,i) = B(:,i)/st; %B(:,i)/((i-1)*sr*2/dy)^(i-1);
    st = st*sr*2/dy*(1+1/(i-1))^(i-1)*i;
end
