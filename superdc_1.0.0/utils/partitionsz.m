function m = partitionsz(D,tr)

m = [];
ch = child(tr);
k = length(tr);
for i = 1:k
    if isempty(ch{i})
        m = [m, size(D{i},1)];
    end
end

