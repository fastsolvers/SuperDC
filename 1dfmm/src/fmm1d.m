%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  fmm w/ adaptive domain partition  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [z, nflops] = fmm1d( r, x, y, q, fun, scaling)

%%%%%%%%%%% Input %%%%%%%%%%%%

% r: truncation order, positive integer
% x: target points, real vectors
% y: source point, real vectors
% q: charges, vectors
% fun:  switch fun
%             case 1
%                 D = 1 ./ (x - y.');
%             case 2
%                 D = 1 ./ (x - y.').^2;
%             case 3
%                 D = log(abs((x - y.')));
%         end


%%%%%%%%%%% Output %%%%%%%%%%%%

% z: Ker(x,y)*q




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% set up %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert x, y to column vectors
if size( x,1 ) < size( x,2 )
    x = x.';
end
if size( y,1 ) < size( y,2 )
    y = y.';
end

% sort x, y
x = sort( x );
y = sort( y );

% partition threshold, 2^9 is best here
N0 = 2^9;

% separation ratio
tau = 0.6;

% stability scaling
if nargin  < 6
    scaling = 1;
end

% computation cell
z2 = max([x;y]);
z1 = min([x;y]);
z2 = z2 + 0.1*abs(z2);
z1 = z1 - 0.1*abs(z1);

nflops = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  adaptive partition %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [center, radius, left_end_point, right_end_point, tree_index, parent, level]
I = [(z1+z2)/2, (z2-z1)/2, 1, 1, 1, 0, 1];
S = {};
S = push(S, I);
Px = {[1:length(x)]};
Py = {[1:length(y)]};

while ~isempty(S)
    [S, i] = pop(S);
    c = i(1); d = i(2); el = i(3); er = i(4); idx = i(5); lvl = i(7);
    n = size(I, 1);
    px = Px{idx}; Px{idx} = [];
    py = Py{idx}; Py{idx} = [];
    
    %%% bisect the interval
    i1 = [c-d/2, d/2, el, 1, n+1, idx, lvl+1];
    i2 = [c+d/2, d/2, 0, er, n+2, idx, lvl+1];
    I = [I; i1; i2];
    if el == 1
        Px{n+1} = px(c-d <= x(px) & x(px) <= c);
        Py{n+1} = py(c-d <= y(py) & y(py) <= c);
    elseif el == 0
        Px{n+1} = px(c-d < x(px) & x(px) <= c);
        Py{n+1} = py(c-d < y(py) & y(py) <= c);
    end
    
    if er == 1
        Px{n+2} = px(c < x(px) & x(px) <= c+d);
        Py{n+2} = py(c < y(py) & y(py) <= c+d);
    elseif er == 0
        Px{n+2} = px(c < x(px) & x(px) < c+d);
        Py{n+2} = py(c < y(py) & y(py) < c+d);
    end
    
    if length(Px{n+1}) > N0 || length(Py{n+1}) > N0
        S = push(S, i1);
    end
    
    if length(Px{n+2}) > N0 || length(Py{n+2}) > N0
        S = push(S, i2);
    end 
end


% tree information
tr = I(:, 6);
ch = child1(tr);
n = length(tr);
post_ord = postorder(1, tr, ch);
pre_ord = preorder(1, tr, ch);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% neighbors and interaction list %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

neighbor = cell(n,1);
interlist = cell(n,1);
for i = pre_ord
    if tr(i) == 0
        neighbor{i} = [];
        interlist{i} = [];
    else
        neighbor{i} = [sib(tr,ch,i)];
        
        cousins = [];
        for j = neighbor{tr(i)}
            if isempty(ch{j})
                cousins = [cousins, j];
            else
                cousins = [cousins, ch{j}(1), ch{j}(2)];
            end
        end
        
        ci = I(i, 1); 
        ri = I(i, 2);
        for j = cousins
            cj = I(j, 1); 
            rj = I(j, 2);
            if ri + rj <= tau * abs(ci - cj)
                interlist{i} = [interlist{i}, j];
            else
                neighbor{i} = [neighbor{i}, j];
            end
        end
    end
end


for i = pre_ord
    if isempty(ch{i})
        for j = neighbor{i}
            neighbor{j} = unique([neighbor{j}, i]);
        end
    end 
    
    for j = interlist{i}
        interlist{j} = unique([interlist{j}, i]);
    end 
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% fmm product  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = cell(n, 1);    % local expansion
v = cell(n, 1);    % multipole expansion
z = zeros(size(x,1), size(q,2));    % K * q


%%%% post order (bottom up) traversal
for i = post_ord
    lvl = I(i, 7);
    if lvl >= 2
        if isempty(ch{i})
            yi = y(Py{i});
            Vi = computeU_scaled(yi, r, I(i,1), 2*I(i,2), scaling);
            v{i} = Vi' * q(Py{i}, :);
            nflops = nflops + size(Vi,2)*(2*size(Vi,1)-1)*size(q,2);
        else
            c1 = ch{i}(1); 
            c2 = ch{i}(2);
            Wc1 = computeT_scaled( r, I(c1,1), I(i,1), 2*I(c1,2), 2*I(i,2), scaling);
            Wc2 = computeT_scaled( r, I(c2,1), I(i,1), 2*I(c2,2), 2*I(i,2), scaling);
            v{i} = Wc1' * v{c1} + Wc2' * v{c2};
            nflops = nflops + flops('prod', Wc1, 't', v{c1}, 'n')...
                                         + flops('prod', Wc2, 't', v{c2}, 'n') + numel(v{i});
        end
    end
end

%%%% pre order (top down) traversal
for i = pre_ord
    lvl = I(i, 7);
    if lvl >= 2
        list = interlist{i};
        for j = list
            Bij = computeB_scaled(r, I(i,1), I(j,1), 2*I(i,2), 2*I(j,2), fun, scaling);
            if isempty(u{i})
                u{i} = Bij * v{j};
                nflops = nflops + flops('prod', Bij, 'n', v{j}, 'n');
            else
                u{i} = u{i} + Bij * v{j};
                nflops = nflops + flops('prod', Bij, 'n', v{j}, 'n') + numel(u{i});
            end
        end
        
        if lvl > 3
            p = tr(i);
            if ~isempty(u{p})
                Ri = computeT_scaled(r, I(i,1), I(p,1),  2*I(i,2), 2*I(p,2), scaling);
                if isempty(u{i})
                    u{i} = Ri * u{p};
                    nflops = nflops + flops('prod', Ri, 'n', u{p}, 'n');
                else
                    u{i} = u{i} + Ri * u{p};
                    nflops = nflops + flops('prod', Ri, 'n', u{p}, 'n') + numel(u{i});
                end
            end
        end
    end
end


%%%%% near field interaction & final product
for i = post_ord
    if isempty(ch{i})
        px = Px{i}; 
        if ~isempty(px)
            xi = x(px);
            %%% evaluate local expansion
            if ~isempty(u{i})
                Ui = computeU_scaled(xi, r, I(i,1), 2*I(i,2), scaling);
                z(px, :) = Ui * u{i};
                nflops = nflops + flops('prod', Ui, 'n', u{i}, 'n');
            end


            %%% near field evaluation
            for j = union(i, neighbor{i})
                py = Py{j};
                if ~isempty(py)
                    yj = y(py);
                    switch fun
                        case 1
                            D = 1 ./ (xi - yj.');
                        case 2
                            D = 1 ./ (xi - yj.').^2;
                        case 3
                            D = log(abs((xi - yj.')));
                    end
                    D(isinf(D)) = 0;
                    z(px, :) = z(px, :) + D * q(py, :);
                    nflops = nflops + size(D,1)*(2*size(D,2)-1)*size(q,2) + numel(px)*size(z,2);
                end
            end
            
        end
    end
end



