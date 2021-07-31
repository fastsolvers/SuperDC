%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% symmetric HSS matrix vector product routine  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, nflops] = symhssmatvec2(D, U, R, B, tr, x)

nflops = double(0);
if nargin < 9 
    t = 0;
end
n = length(tr); 
u = cell(n,1); 
v = cell(n,1); 
ch = child(tr);
m = [];
for i = 1:n
    if ~(i > length(U))
        if ~isempty(U{i})
            m = [m, size(U{i},1)];
        end
    end
end
rg = indrange(tr,m);


if n > 1
    y = zeros(sum(m), size(x,2));
else
    y = D{1} * x;
    nflops = nflops + flops('prod', D{1}, 'n', x, 'n');
    return
end

%%% compute mat-vec without transpose
% bottom up traversal (multipole expansion)
for i = 1:n-1
    if isempty(ch{i})
        v{i} = U{i}' * x(rg(i,1):rg(i,2), :);
        nflops = nflops + flops('prod', U{i}, 't', x(rg(i,1):rg(i,2), :), 'n');
        
        if ~isempty(D{i})
            y(rg(i,1):rg(i,2), :) = D{i} * x(rg(i,1):rg(i,2), :);
            nflops = nflops + flops('prod', D{i}, 'n', x(rg(i,1):rg(i,2), :), 'n');
        end
    else
        v{i} = R{ch{i}(1)}' * v{ch{i}(1)} + R{ch{i}(2)}' * v{ch{i}(2)};
        nflops = nflops + flops('prod', R{ch{i}(1)}, 't', v{ch{i}(1)}, 'n')...      
                                      + flops('prod', R{ch{i}(2)}, 't', v{ch{i}(2)}, 'n') + numel(v{i});
    end
end

% top down traversal (local expansion)
for i = n-1:-1:1
    if isempty(ch{i})
        if i == ch{tr(i)}(1)
            if tr(i) ~= n
                y(rg(i,1):rg(i,2), :) = y(rg(i,1):rg(i,2), :) + U{i} * (R{i} * u{tr(i)} + B{i} * v{ch{tr(i)}(2)});
                nflops = nflops + flops('prod', R{i}, 'n', u{tr(i)}, 'n') + flops('prod', B{i}, 'n', v{ch{tr(i)}(2)}, 'n')...
                                              + size(R{i},1) * size(u{tr(i)},2) + size(U{i},1) * (2*size(U{i},2)-1) * size(u{tr(i)},2);
            else
                y(rg(i,1):rg(i,2), :) = y(rg(i,1):rg(i,2), :) + U{i} * B{i} * v{ch{tr(i)}(2)};
                nflops = nflops + flops('prod', B{i}, 'n', v{ch{tr(i)}(2)}, 'n')...
                                              + size(U{i},1) * (2*size(U{i},2)-1) * size(v{ch{tr(i)}(2)},2);
            end
        else
            sib = ch{tr(i)}(1);
            if tr(i) ~= n
                y(rg(i,1):rg(i,2), :) = y(rg(i,1):rg(i,2), :) + U{i} * ( R{i} * u{tr(i)} + B{sib}' * v{ch{tr(i)}(1)} );
                nflops = nflops + flops('prod', R{i}, 'n', u{tr(i)}, 'n') + flops('prod', B{sib}, 't', v{ch{tr(i)}(1)}, 'n')...
                                              + size(R{i},1) * size(u{tr(i)},2) + size(U{i},1) * (2*size(U{i},2)-1) * size(u{tr(i)},2);
            else
                y(rg(i,1):rg(i,2), :) = y(rg(i,1):rg(i,2), :) + U{i} * B{sib}' * v{ch{tr(i)}(1)};
                nflops = nflops + flops('prod', B{sib}, 't', v{ch{tr(i)}(1)}, 'n')...
                                              + size(U{i},1) * (2*size(U{i},2)-1) * size(v{ch{tr(i)}(1)},2);
            end
        end
    else
        if i == ch{tr(i)}(1)
           if tr(i) == n
                u{i} = B{i} * v{ch{tr(i)}(2)};
                nflops = nflops + flops('prod', B{i}, 'n', v{ch{tr(i)}(2)}, 'n');
           else 
                u{i} = R{i} * u{tr(i)} + B{i} * v{ch{tr(i)}(2)};
                nflops = nflops +  flops('prod', R{i}, 'n', u{tr(i)}, 'n')...
                                              + flops('prod', B{i}, 'n', v{ch{tr(i)}(2)}, 'n') + numel(u{i});
           end
        else
            sib = ch{tr(i)}(1);
            if tr(i) == n
                u{i} = B{sib}' * v{ch{tr(i)}(1)};
                nflops = nflops + flops('prod', B{sib}, 't', v{ch{tr(i)}(1)}, 'n');
            else 
                u{i} = R{i} * u{tr(i)} + B{sib}' * v{ch{tr(i)}(1)};
                nflops = nflops +  flops('prod', R{i},'n', u{tr(i)}, 'n')...
                                              + flops('prod', B{sib}, 't', v{ch{tr(i)}(1)}, 'n') + numel(u{i});
            end
        end
    end
end



