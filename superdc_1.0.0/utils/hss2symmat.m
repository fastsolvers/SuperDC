function A = hss2symmat(D, U, R, B, tr)


m = partitionsz(D,tr);
N = sum(m);
A = zeros(N);
ch = child(tr);
rg = indrange(tr, m);
n = length(tr);

SU = {};
for i = 1:n
    if isempty(ch{i})
        SU = push(SU, U{i});
        A(rg(i,1):rg(i,2), rg(i,1):rg(i,2)) = D{i};
        continue;
    end
    
    c1 = ch{i}(1);
    c2 = ch{i}(2);
    
    [SU, Uc2] = pop(SU);
    [SU, Uc1] = pop(SU);
    
    A(rg(c1,1):rg(c1,2), rg(c2,1):rg(c2,2)) = Uc1 * B{c1} * Uc2';
    A(rg(c2,1):rg(c2,2), rg(c1,1):rg(c1,2)) = A(rg(c1,1):rg(c1,2), rg(c2,1):rg(c2,2))';

    if i < n
        Ui = [Uc1 * R{c1}; Uc2 * R{c2}];
        SU = push(SU, Ui);
    end
end