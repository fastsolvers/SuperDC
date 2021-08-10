function [D,U,R,B,W,V] = band2hss(A,w,tr,m)
% HSS form of a banded matrix A
% Input:
% A:  banded form
% w:  half bandwidth, total bandwidth = 2*w+1
% tr: postordered binary tree, e.g., tr = [3 3 7 6 6 7 0];
% m:  partition, e.g., m = [100 100 100 100] when n = 400
% Output:
% D,U,etc.: HSS of the banded form
% Last updated: 2016/11/5.  J. Xia
% Last updated: 2016/11/11. J. Xia
% Last updated: 2016/12/4. J. Xia

D = {[]}; U = {[]}; R = {[]}; B = {[]}; W = {[]}; V = {[]};

n = length(tr);

if n == 1
    D{1} = T;
    return;
end

ch = child(tr);

lm = length(m);

bn = 1:n;

cl = [];
for i = 1:n
    if isempty(ch{i})
        cl = [cl i];
    end
end

l(1) = 1;
for i = 1:lm
    l(i+1) = l(i)+m(i);
end
O = zeros(w,w); I = eye(w);

k = 0;
for i = 1:n-1
    if isempty(ch{i})
        k = k+1;
        D{i} = A(l(k):l(k+1)-1,l(k):l(k+1)-1);
        U{i} = zeros(m(k),2*w); U{i}(1:w,1:w) = I; U{i}(end-w+1:end,w+1:2*w) = I;
        V{i} = U{i};
    end
    
    if ch{tr(i)}(1) == i
        rc = i;
        while ~isempty(ch{rc})
            rc = ch{rc}(2);
        end
        ii = find(cl == rc);
        B{i} = [O O; A(l(ii+1)-w:l(ii+1)-1,l(ii+1):l(ii+1)+w-1) O];
        R{i} = [I O;O O];
    else
        rc = ch{tr(i)}(1);
        while ~isempty(ch{rc})
            rc = ch{rc}(2);
        end
        ii = find(cl == rc);
        B{i} = [O A(l(ii+1):l(ii+1)+w-1,l(ii+1)-w:l(ii+1)-1); O O];
        R{i} = [O O;O I];
    end
    W{i} = R{i};
end

function ch = child(tr)

ch = {[]};

for i = 1:length(tr)-1
    ch{tr(i)} = [];
end

for i = 1:length(tr)-1
    ch{tr(i)} = [ch{tr(i)} i];
end