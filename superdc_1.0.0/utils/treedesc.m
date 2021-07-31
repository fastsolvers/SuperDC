function [desc1, desc2] = treedesc(tr)
% desc1: smallest leaf descendant of each node
% desc2: largest leaf descendant of each node

n = length(tr);
ch = child(tr);

desc1 = zeros(1,n);
desc2 = zeros(1,n);
for i = 1:n
    if isempty(ch{i})
        desc1(i) = i;
        desc2(i) = i;
    else
        desc1(i) = desc1(ch{i}(1));
        desc2(i) = desc2(ch{i}(2));
    end
end