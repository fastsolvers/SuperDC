   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%  convert a block-tri-diagonal matrix into hss form  %%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [D, U, R, B, W, V, tr, m] = blktri2hss(D1, F1, G1, N)
%%% Input:
%%% D1: cell array of diagonal blocks
%%% F1: cell array of (-1) subdiagonal blocks
%%% G1: cell array of (+1) superdiagonal blocks
%%% N: number of blocks that constitute the hss leaf block

%%% Output:
%%% D, U, R, B, W, V: hss generators
%%% tr: hss tree
%%% m: hss leaf blocks size

if nargin < 4
    N = 3;
end

% size of each small blocks
num_blks = length(D1);
m1 = zeros(num_blks, 1);
for j = 1:length(m1)
    m1(j) = size(D1{j}, 1);
end


J = 1:N:num_blks;
K = N * ones(1, length(J)); 
K(end) = num_blks - J(end) + 1;


% size of each leaf node
m = zeros(length(J), 1);
for j = 1:length(J)
    for i = 0:K(j)-1
        m(j) = m(j) + size(D1{J(j)+i}, 1);
    end
end
if m(end) <= m(end-1)/2
    m(end-1) = m(end-1) + m(end);
    m = m(1:end-1);
    J = J(1:end-1);
    K(end-1) = K(end-1) + K(end);
    K = K(1:end-1);
end

%%% off diagonal blocks of each hss leaf nodes
F = cell(length(m), 1);
G = cell(length(m), 1);
for j = 1:length(m)
    if j == 1
        G{j} = G1{J(j) + N - 1};
    elseif j == length(m)
        F{j} = F1{J(j)};
    else
        G{j} = G1{J(j) + N - 1};
        F{j} = F1{J(j)};
    end
end


%%% hss tree info
tr = n2tree(2 * length(J) - 1);
ch = child(tr);
[desc1, desc2] = treedesc(tr);
k = length(tr);
[~, ~, ~, leaf] = lvlnodes(tr);

%%% initialization
D = cell(k,1); U = cell(k,1); R = cell(k,1);
B = cell(k,1); W = cell(k,1); V = cell(k,1);


j = 0;       % leaf counter
for i = 1:k
    if isempty(ch{i})
        j = j + 1;
        %%% D generator at leaf node
        D{i} = zeros(m(j));
        idx = 1;
        for p = 0:K(j)-1
            curr_blk = J(j)+p;
            DD = D1{curr_blk};
            D{i}(idx:idx+m1(curr_blk)-1, idx:idx+m1(curr_blk)-1) = DD;
            if p ~= 0
                A1 = F1{curr_blk};
                D{i}(idx:idx+size(A1,1)-1, idx-size(A1,2):idx-1) = A1;
            end
            idx = idx + m1(curr_blk);
            if p ~= K(j)-1
                A2 = G1{curr_blk};
                D{i}(idx-size(A2,1):idx-1, idx:idx+size(A2,2)-1) = A2;
            end   
        end
        
        %%% U, V generator at leaf nodes
        if j == 1
            U{i} = [ zeros(m(j) - size(G{j},1), size(G{j},1));
                                                              eye(size(G{j},1)) ];
            V{i} = [ zeros(m(j) - size(F{j+1},2), size(F{j+1},2));
                                                              eye(size(F{j+1},2)) ];
        elseif j == (k+1)/2
            U{i} = [                                    eye(size(F{j},1));  
                         zeros(m(j) - size(F{j},1), size(F{j},1)) ];
            V{i} = [                                    eye(size(G{j-1},2));  
                         zeros(m(j) - size(G{j-1},2), size(G{j-1},2)) ];
        else
            U{i} = [ eye(size(F{j},1)),      zeros(size(F{j},1), size(G{j},1));
                          zeros(m(j)-size(F{j},1)-size(G{j},1), size(F{j},1)+size(G{j},1));
                          zeros(size(G{j},1), size(F{j},1)),     eye(size(G{j},1)) ];
                     
            V{i} = [ eye(size(G{j-1},2)),      zeros(size(G{j-1},2), size(F{j+1},2));
                          zeros(m(j)-size(G{j-1},2)-size(F{j+1},2), size(G{j-1},2)+size(F{j+1},2));
                          zeros(size(F{j+1},2), size(G{j-1},2)),     eye(size(F{j+1},2)) ];
        end
        
        
        
        
    else     
        c1 = ch{i}(1); c2 = ch{i}(2);
        js_c1 = desc1(c1); js_c2 = desc1(c2);
        jl_c1 = desc2(c1); jl_c2 = desc2(c2);
        js_c1_idx = find(leaf == js_c1);
        js_c2_idx = find(leaf == js_c2);
        jl_c1_idx = find(leaf == jl_c1);
        jl_c2_idx = find(leaf == jl_c2);
        
        %%% R, B, W generators
        if js_c1_idx == 1 && jl_c2_idx ~= (k+1)/2
            R{c1} = [ zeros(size(G{jl_c1_idx}, 1), size(G{jl_c2_idx}, 1)) ];
            
            R{c2} = [ zeros(size(F{js_c2_idx}, 1), size(G{jl_c2_idx}, 1));
                                                                     eye(size(G{jl_c2_idx}, 1)) ];
                                                                 
            W{c1} = [ zeros(size(F{js_c2_idx}, 2), size(F{jl_c2_idx+1}, 2)) ];
            
            W{c2} = [ zeros(size(G{jl_c1_idx}, 2), size(F{jl_c2_idx+1}, 2));
                                                                            eye(size(F{jl_c2_idx+1}, 2)) ];
            
            B{c1} = [ G{jl_c1_idx},  zeros(size(G{jl_c1_idx}, 1), size(F{jl_c2_idx+1}, 2)) ];
            
            B{c2} = [ F{js_c2_idx};  
                            zeros(size(G{jl_c2_idx }, 1), size(F{js_c2_idx}, 2)) ];
            
        elseif js_c1_idx ~=1 && jl_c2_idx == (k+1)/2
            R{c1} = [                                            eye(size(F{js_c1_idx}, 1));
                            zeros(size(G{jl_c1_idx}, 1), size(F{js_c1_idx}, 1)) ];
                        
            R{c2} = [ zeros(size(F{js_c2_idx}, 1), size(F{js_c1_idx}, 1)) ];
            
            W{c1} = [                                          eye(size(G{js_c1_idx-1}, 2));
                             zeros(size(F{js_c2_idx}, 2), size(G{js_c1_idx-1}, 2)) ];
                        
            W{c2} = [ zeros(size(G{jl_c1_idx}, 2), size(G{js_c1_idx-1}, 2)) ];
            
            B{c1} = [ zeros(size(F{js_c1_idx}, 1), size(G{jl_c1_idx}, 2));   
                            G{jl_c1_idx} ];
                        
            B{c2} = [ zeros(size(F{js_c2_idx}, 1), size(G{js_c1_idx-1}, 2)),  F{js_c2_idx} ];
            
        elseif js_c1_idx ~= 1 && jl_c2_idx ~= (k+1)/2
            R{c1} = [ eye(size(F{js_c1_idx}, 1)), zeros(size(F{js_c1_idx}, 1), size(G{jl_c2_idx}, 1));
                             zeros(size(G{jl_c1_idx}, 1), size(F{js_c1_idx}, 1)+size(G{jl_c2_idx}, 1)) ];
                        
            R{c2} = [ zeros(size(F{js_c2_idx}, 1), size(F{js_c1_idx}, 1)+size(G{jl_c2_idx}, 1));
                             zeros(size(G{jl_c2_idx}, 1), size(F{js_c1_idx}, 1)), eye(size(G{jl_c2_idx}, 1)) ];
                                                                 
            W{c1} = [ eye(size(G{js_c1_idx-1}, 2)), zeros(size(G{js_c1_idx-1}, 2), size(F{jl_c2_idx+1}, 2));
                              zeros(size(F{js_c2_idx}, 2), size(G{js_c1_idx-1}, 2)+size(F{jl_c2_idx+1}, 2)) ];
                         
            W{c2} = [ zeros(size(G{jl_c1_idx}, 2), size(G{js_c1_idx-1}, 2)+size(F{jl_c2_idx+1}, 2));
                              zeros(size(F{jl_c2_idx+1}, 2), size(G{js_c1_idx-1}, 2)), eye(size(F{jl_c2_idx+1}, 2)) ];
                                                                        
            B{c1} = [ zeros(size(F{js_c1_idx}, 1), size(G{jl_c1_idx}, 2)+size(F{jl_c2_idx+1}, 2));
                             G{jl_c1_idx},                 zeros(size(G{jl_c1_idx},1), size(F{jl_c2_idx+1}, 2)) ];
                         
            B{c2} = [ zeros(size(F{js_c2_idx}, 1), size(G{js_c1_idx-1}, 2)),                F{js_c2_idx};
                             zeros(size(G{jl_c2_idx}, 1), size(G{js_c1_idx-1}, 2)+size(F{js_c2_idx}, 2)) ];
            
        elseif js_c1_idx == 1 && jl_c2_idx == (k+1)/2
            B{c1} = G{jl_c1_idx};
            B{c2} = F{js_c2_idx};
            
        end  
    end
end