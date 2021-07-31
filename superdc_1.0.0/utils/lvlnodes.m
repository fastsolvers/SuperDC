function [node, leftnode, rightnode, leaf] = lvlnodes(tr)

ch = child(tr);
lvl = hsslevel(tr);
k = max(lvl);
n = length(tr);

node = cell(k,1); 
rightnode = cell(k,1); 
leftnode = cell(k,1);
leaf = [];

%%% left and right nodes of each level
for i = 1:n-1
    node{lvl(i)} = [node{lvl(i)}, i];
    
    if isempty(ch{i})
        leaf = [leaf, i];
    end

    if i == ch{tr(i)}(1)
        leftnode{lvl(i)} = [leftnode{lvl(i)}, i];
    else
        rightnode{lvl(i)} = [rightnode{lvl(i)}, i];
    end
end