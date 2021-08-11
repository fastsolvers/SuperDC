%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SuperDC dividing stage (old, unstable) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [D, B, Z, desc, nflops] = divide1(D, U, B, R, tr)
%%% Input:
%%% D, U, B, R: hss generators
%%% tr: hss tree

%%% Output:
%%% D: D generators after dividing
%%% B: B generators after dividing
%%% Z: rank-r update at each nonleaf node

n = length(tr);
ch = child(tr);
nflops = 0;

% smallest descendant of each node
desc = treedesc(tr);

%%% dividing process
for i = n:-1:1
    if isempty(ch{i})
        continue;
    end

    
    % divide at non-leaf nodes
    c1 = ch{i}(1); 
    c2 = ch{i}(2);


    % update B generators of left subtree
    if isempty(ch{c1})
        if size(B{c1}, 1) <= size(B{c1}, 2)
            D{c1} = D{c1} - (U{c1} * U{c1}');
            nflops = nflops + flops('prod', U{c1}, 'n', U{c1}, 't') + 2 * flops('sum', D{c1});
        else
            T = U{c1} * B{c1};
            nflops = nflops + flops('prod', U{c1}, 'n', B{c1}, 'n');
            D{c1} = D{c1} - (T * T');
            nflops = nflops + flops('prod', T, 'n', T, 't') + 2 * flops('sum', D{c1});
        end
    else
        if size(B{c1}, 1) <= size(B{c1}, 2)
            Sp = eye(size(R{ch{c1}(1)}, 2)); 
        else
            Sp = B{c1};
        end
        S = {};
        S = push(S, Sp);                   
            
        for j = c1-1 : -1 : desc(c1)
            [S, Sp] = pop(S);
            if j == ch{tr(j)}(2)            
                sib = ch{tr(j)}(1);
                S = push(S, Sp);       
                B{sib} = B{sib} - R{sib} * (Sp * Sp') * R{j}';
                nflops = nflops + flops('prod', Sp, 'n', Sp, 't') + flops('prod', R{j}, 'n', Sp*Sp', 'n')...
                    + flops('prod', R{j} * (Sp * Sp'), 'n', R{sib}, 't') + flops('sum', B{sib});
            end
          
            Sj = R{j} * Sp;
            nflops = nflops + flops('prod', R{j}, 'n', Sp, 'n');
            if isempty(ch{j})
                T = U{j} * Sj;
                nflops = nflops + flops('prod', U{j}, 'n', Sj, 'n');
                D{j}  = D{j} - T * T';
                nflops = nflops + flops('prod', T, 'n', T, 't') + flops('sum', D{j});
            else 
                S = push(S, Sj);
            end    
        end
    end
    
    
    % update B generators of right subtree
    if isempty(ch{c2})
        if size(B{c1}, 1) <= size(B{c1}, 2)
            T = U{c2} * B{c1}';
            nflops = nflops + flops('prod', U{c2}, 'n', B{c1}, 'T');
            D{c2} = D{c2} - (T * T');
            nflops = nflops + flops('prod', T, 'n', T, 't') + 2 * flops('sum', D{c2});
        else
            D{c2} = D{c2} - (U{c2} * U{c2}');
            nflops = nflops + flops('prod', U{c2}, 'n', U{c2}, 't') + 2 * flops('sum', D{c2});
        end
    else
        if size(B{c1}, 1) <= size(B{c1}, 2)
            Sp = B{c1}';
        else
            Sp = eye(size(R{ch{c2}(1)}, 2)); 
        end
        S = {};
        S = push(S, Sp);
        
        for j = c2-1 : -1 : desc(c2)
            [S, Sp] = pop(S);
            if j == ch{tr(j)}(2)   
                S = push(S, Sp);               
                sib = ch{tr(j)}(1);
                B{sib} = B{sib} - R{sib} * (Sp * Sp') * R{j}';
                nflops = nflops + flops('prod', Sp, 'n', Sp, 't') + flops('prod', R{j}, 'n', Sp*Sp', 'n')...
                    + flops('prod', R{j} * (Sp * Sp'), 'n', R{sib}, 't') + flops('sum', B{sib});  
            end
            
            Sj = R{j} * Sp;
            nflops = nflops + flops('prod', R{j}, 'n', Sp, 'n');
             if isempty(ch{j})
                 T = U{j} * Sj;
                 nflops = nflops + flops('prod', U{j}, 'n', Sj, 'n');
                 D{j}  = D{j} - T * T';
                 nflops = nflops + flops('prod', T, 'n', T, 't') + flops('sum', D{j});
             else
                 S = push(S, Sj);
             end 
        end
    end
end


%%% rank-r updates at non-leaf nodes
Z = cell(n,1);
SU = {};
for i = 1:n
    if isempty(ch{i})
        SU = push(SU, U{i});
        continue;
    end
    
    c1 = ch{i}(1);
    c2 = ch{i}(2);
    
    [SU, Uc2] = pop(SU);
    [SU, Uc1] = pop(SU);
    if size(B{c1}, 1) <= size(B{c1}, 2)
        Z{i} = [Uc1; 
                    Uc2 * B{c1}'];
        nflops = nflops + flops('prod', Uc2, 'n', B{c1}, 't');
    else
        Z{i} = [Uc1 * B{c1};
                    Uc2];
        nflops = nflops + flops('prod', Uc1, 'n', B{c1}, 'n');
    end
    
    if i < n
        Ui = [Uc1 * R{c1}; Uc2 * R{c2}];
        nflops = nflops + flops('prod', Uc1, 'n', R{c1}, 'n') + flops('prod', Uc2, 'n', R{c2}, 'n');
        SU = push(SU, Ui);
    end
end
