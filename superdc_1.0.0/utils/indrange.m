function rg = indrange(tr,m)
%%%%% index range of each node 

n = length(tr);
rg = zeros(n,2); 
ch = child(tr);

rg(1,:) = [1 m(1)];
lt = 1; it = 1;
if n > 1
    for i = 1:n
        if isempty(ch{i})
            rg(i,:) = [lt lt+m(it)-1];
            lt = rg(i,2)+1; it = it+1;
        else
            rg(i,:) = [rg(ch{i}(1),1) rg(ch{i}(2),2)];
        end
    end
end