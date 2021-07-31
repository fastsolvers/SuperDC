function storage = Qstorage(Q)
%%% Q: HSS eigenmatrix
%%% tr: HSS tree

[Q0, I, tr, m] = Q{:};
ch = child(tr);
k = length(tr);
storage = numel(I) + numel(m) + numel(tr);

for i = 1:k
    if isempty(ch{i})
        storage = storage + numel(Q0{i});
    else
        r = length(Q0{i});
        for j = 1:r
            for l = 1:6
                storage = storage + numel(Q0{i}{j}{1}{l});
            end
            
            for l = 2:7
                storage = storage + numel(Q0{i}{j}{l});
            end
        end
    end
    
end