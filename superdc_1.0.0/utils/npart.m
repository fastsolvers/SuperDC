function [tr,m] = npart(n,ni,tr,ltr)

if nargin == 2 || ltr == 0
    k = floor(n/ni);
    m = ni*ones(1,k);
    mr = mod(n,ni);
    if mr >= ni/2
        k = k+1;
        m = [m mr];
    else
        m(end) = m(end)+mr;
    end
    tr = n2tree(2*k-1);
else
    k = (ltr+1)/2;
    ni = floor(n/k);
    m = ni*ones(1,k);
    m(end) = m(end)+mod(n,ni);
end