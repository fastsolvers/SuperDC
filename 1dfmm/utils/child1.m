function ch = child1(tr)

n = length(tr);
ch = cell(1,n);

for i = 2:n
    ch{tr(i)} = [];
end

for i = 2:n
    ch{tr(i)} = [ch{tr(i)} i];
end