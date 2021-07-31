%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% root of the secular equation %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function  [x, tau,org, nflops, percent] = rootfinder(d, v, N)
%%% Input:
%%% d, v: as in secular equation
%%% N: size threshold to use fmm

%%% Ouput:
%%% x: roots of secular equation 


if nargin < 3
    N = 2^13;
end
nflops = 0;
n = length(v);
d = reshape(d, [n, 1]);
v = reshape(v, [n, 1]);
if n == 1
    x = v^2 + d;
    tau = x - d;
    org = 1;
    nflops = 2;
    return
end

C = 64;
displaywarning = 1;
record = 1;
stop_criteria = 1;
alpha = 0;
r = 50;
N0 = 2^20;
percent = 1;

FMM_ITER = max(ceil(log(n)/log(2))-6, 5);
FLAG = 1;
MAX_ITER = 100;

if record
    iter = zeros(n,1);
    residual = zeros(n,1);
    erretms = zeros(n,1);
end

v_norm = norm(v);
rho = 1 / v_norm^2;
v = v / norm(v);
d0 = diff(d);
x = (d(1:n-1) + d(2:n)) / 2;
org = [1:n-1].';

if n >= N
    [fl, fu, nflops1] = trifmm1d_local_shift(r, x, d, v.^2, d0/2, org, 1);
    f0 = rho - fl - fu;
%     [f0, nflops1] = fmm1d_local_shift(r, x, d, v.^2, d0/2, org, 1);
%     f0 = rho - f0;
    nflops = nflops + nflops1;
else
    K = d(org) - d.';
    K = bsxfun(@plus, K, d0/2);
    K = 1 ./ K;
    f0 = rho - K * v.^2;
    nflops = nflops + flops('prod', K, 'n', v, 'n');
end
h = 2 * diff(v.^2) ./ d0;
g = f0 - h;
org = [1:n-1]' + (f0<0);


%%% upper bound and lower bound 
d1 = d(1:n-1); 
d2 = d(2:n);
v1 = v(1:n-1);
v2 = v(2:n);
xub = zeros(n-1,1); 
xlb = zeros(n-1,1);
I1 = find(f0>=0);
I2 = find(f0<0);
xub(I1) = ( d2(I1) - d1(I1) ) / 2; 
xlb(I2) = ( d1(I2) - d2(I2) ) / 2;
xlb0 = xlb;
xub0 = xub;


%%% initial guess
a = (2*(f0>=0) - 1) .* d0 .* g + v(1:n-1).^2 + v(2:n).^2;
b = ((f0>=0) .* v(1:n-1).^2 - (f0<0) .* v(2:n).^2) .* d0;
tau = zeros(n-1,1);
I1 = find(a>0); 
I2 = find(a<=0);
tau(I1) = 2*b(I1) ./ (a(I1) + sqrt(abs(a(I1).^2 - 4*b(I1).*g(I1))));
tau(I2) = (a(I2) - sqrt(abs(a(I2).^2 - 4*b(I2).*g(I2)))) ./ (2*g(I2));


%%% bound check
I = find((tau<=xlb) | (tau>=xub));
if ~isempty(I)
    tau(I) = (xub(I) + xlb(I)) / 2;
end
x = tau + d(org);
nflops = nflops + 18*n;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  iterate n-1 roots simultaneously  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iter_ct = 0;
swtch = ones(n-1,1);
while (iter_ct < FMM_ITER) && FLAG
    if iter_ct > 0
        pref = f;
    end

    
    %%% bound check with bisection safeguard
    I = find(tau==0 | tau.*xub<0 | tau.*xlb<0 | isnan(tau) | tau<xlb0 | tau>xub0);
    if ~isempty(I)
        tau(I) = (xub(I) + xlb(I)) / 2;
        x(I) = tau(I) + d(org(I));
    end
    
    
    %%% secular function evaluation
     if n >= N
        [psi, phi, nflops1] = trifmm1d_local_shift(r, x, d, v.^2, tau, org, 1);
        f = rho - psi - phi;
        nflops = nflops + nflops1;
        
        [dpsi, dphi, nflops1] = trifmm1d_local_shift(r, x, d, v.^2, tau, org, 2);
        df = dpsi + dphi;
        nflops = nflops + nflops1;
        
     else   
        K = d(org) - d.';
        K = bsxfun(@plus, K, tau);
        
        K = 1./ K;
        K1 = tril(K);
        K2 = triu(K,1);
        psi = K1 * v.^2;
        phi = K2 * v.^2;
        f = rho - psi - phi;
        nflops = nflops + flops('prod', K, 'n', v, 'n');

        dK = K.^2;
        dK1 = tril(dK);
        dK2 = triu(dK,1);
        dpsi = dK1 * v.^2;
        dphi = dK2 * v.^2;
        df = dpsi + dphi;
        nflops = nflops + flops('prod', dK, 'n', v, 'n');
    end
    
    
    %%% adaptive FMM iterations
     res = rho + abs(psi) + abs(phi);
     nflops = nflops + 2*n;
     if iter_ct > 0
        nJ2 = length(J2);
     end
     J2 = find(abs(f)>C*n*eps*max([res, alpha*ones(n-1,1)], [], 2));
     if iter_ct > 0 && ( length(J2) <= 0.01 * n || abs(nJ2 - length(J2)) < 0.01 * n )
         FLAG = FLAG - 1;
     end
    if iter_ct == 5 || (iter_ct < 5 && FLAG == 0)
        percent = length(J2) / (n-1);
    end
     
    %%% update upper and lower bounds
    I = find(f<0);
    II = find(f>=0);
    xlb(I) = max([tau(I), xlb(I)], [], 2);
    xub(II) = min([tau(II), xub(II)], [], 2);
    
    
    % swtch = 1 : fixed weight method
    % swtch = 0 : middle way method
    if iter_ct > 0
        Iswtch = find((pref .* f > 0) & (abs(pref) > abs(f) / 10));
        swtch(Iswtch) = - swtch(Iswtch) + 1;
    end
    
    
    %%% quadratic equation coefficients.
    a = ((d1 - x) + (d2 - x)) .* f - (d1 - x) .* (d2 - x) .* df;
    b = (d1 - x) .* (d2 - x) .* f;
    if iter_ct > 0
        c = zeros(n-1, 1);
        c0 =  - (d1 - x) .* dpsi - (d2 - x) .* dphi + f;              
        c1 =  - (f0>=0).* ((d2 - x).*df + v1.^2 .* (d1 - d2) ./ (d1 - x).^2)...
                  - (f0<0).* ((d1 - x).*df + v2.^2 .* (d2 - d1) ./  (d2 - x).^2) + f;    
        c(swtch == 0) = c0(swtch == 0);
        c(swtch == 1) = c1(swtch == 1);
    else
        c =  - (f0>=0).* ((d2 - x).*df + v1.^2 .* (d1 - d2) ./ (d1 - x).^2)...
                - (f0<0).* ((d1 - x).*df + v2.^2 .* (d2 - d1) ./  (d2 - x).^2) + f;    
    end
    nflops = nflops + 23*n;
    
     
    %%% eta : update of root
    eta = zeros(n-1,1);
    I1 = find(a>0); 
    I2 = find(a<=0);
    eta(I1) = 2*b(I1) ./ (a(I1) + sqrt(abs(a(I1).^2 - 4*b(I1).*c(I1))));
    eta(I2) = (a(I2) - sqrt(abs(a(I2).^2 - 4*b(I2).*c(I2)))) ./ (2*c(I2));
    Ic = find(abs(c)==0);       % handle corner case : c == 0
    if ~isempty(Ic)
        Ia = find(abs(a(Ic))==0);
        if ~isempty(Ia)
            I = Ic(Ia);
            a(I) = (f0(I)>=0).*(v(I).^2 + (d(I+1)-x(I)).^2 .* (df(I) - v(I).^2 ./ (d(I)-x(I)).^2) )...
                + (f0(I)<0).*(v(I+1).^2 + (d(I)-x(I)).^2 .* (df(I) - v(I+1).^2 ./ (d(I+1)-x(I)).^2) );
        end
        eta(Ic) = b(Ic) ./ a(Ic);
    end
    nflops = nflops + 10*n;


    %%% f*eta should be negative, otherwise run a Newton step
    I = find(f.*eta >= 0);
    if ~isempty(I)
        eta(I) = - f(I) ./ df(I); 
    end
    
    %%% bound check with bisection safeguard
    tmp = tau + eta; 
    I = find( tmp<xlb | tmp>xub | tmp==0 | tmp.*xub<0 | tmp.*xlb<0 | isnan(tmp) );
    if ~isempty(I)
        eta(I) = ( (f(I)>=0).*xlb(I) + (f(I)<0).*xub(I) - tau(I) ) / 2;
    end
    nflops = nflops + numel(I);
    
    %%% update root
    tau = tau + eta;
    x = tau + d(org);
    iter_ct = iter_ct + 1;
end



%%%%%%%%%% check residual after MAX_ITER iterations %%%%%%%%%%
if FMM_ITER
    if n >= N     
        [psi, phi, nflops1] = trifmm1d_local_shift( r, x, d, v.^2, tau, org, 1);
        f = rho - psi - phi;
        nflops = nflops + nflops1;
    else  
        K = d(org) - d.';
        K = bsxfun(@plus, K, tau);
        K = 1 ./ K;
        K1 = tril(K);
        K2 = triu(K,1);
        psi = K1*v.^2;
        phi = K2*v.^2;
        f = rho - psi - phi;
        nflops = nflops + flops('prod',K, 'n', v, 'n');
    end
    residual = abs(f);    
    res = rho + abs(psi) + abs(phi); 
    J2 = find(abs(f)>C*n*eps*max([res, alpha*ones(n-1,1)],[],2));
    nflops = nflops + 9*n;
    if n > N0
        fprintf('rootfinder: continue to iterate on %i roots, out of %i, max resid [%e, %e]\n',...
                                    length(J2), n, max(abs(f(J2))), max(res(J2)))
    end
else
    if n > N0
        fprintf('rootfinder: iterate on %i roots individually\n', n-1)
    end
    J2 = [1:n-1]';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% find roots x(J2) %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% set up %%%
psi = @(x, delta, i) dot(v(1:i).^2, 1./(delta(1:i) - x));
dpsi = @(x, delta, i) dot(v(1:i).^2, 1./(delta(1:i) - x).^2);
phi = @(x, delta, i) dot(v(n:-1:i+1).^2, 1./(delta(n:-1:i+1) - x));
dphi = @(x, delta, i) dot(v(n:-1:i+1).^2, 1./(delta(n:-1:i+1) - x).^2);

QJ = zeros(n,0);
for i = J2'
    y0 = (d(i)+d(i+1))/2;
    w = rho + psi(y0, d, i) + phi(y0, d, i);
    if w >= 0
        org0 = i;                      % origin
        delta = d - d(i);          % shift
        ub = delta(i+1)/2;
        lb = 0;
        if lb < tau(i) && tau(i) <= ub
            x0  = tau(i);
        else
            x0 = (ub + lb) / 2;
        end
        
    elseif w < 0
        org0 = i + 1;                  % origin
        delta = d - d(i+1);        % shift
        ub = 0;
        lb = delta(i)/2;
        if lb < tau(i) && tau(i) < ub
            x0  = tau(i);
        else
            x0 = (ub + lb) / 2;
        end
    end
    
    % initial upper/lower bound
    ub0 = ub;
    lb0 = lb;
    
    %%%% rational interpolation
    psi0 = psi(x0, delta, i); 
    phi0 = phi(x0, delta, i);
    dpsi0 = dpsi(x0, delta, i); 
    dphi0 = dphi(x0, delta, i);
    w = rho + (psi0 + phi0);
    dw = dpsi0 + dphi0;
    nflops = nflops + 16*n;
    if stop_criteria  == 0
        erretm = abs(dot([org0-1:-1:1]', v(1:org0-1).^2 ./ (delta(1:org0-1)-x0))) + ...
                                  dot([n-org0:-1:1]', v(n:-1:org0+1).^2 ./ (delta(n:-1:org0+1)-x0));
        erretm = 8 * (phi(x0, delta, org0) - psi(x0, delta, org0-1)) + erretm + 2 * rho...
                            + 3 * abs(v(org0)^2 / (delta(org0)-x0)) + abs(x0) * dw;
    elseif stop_criteria  == 1
        erretm = C * n * (rho + abs(psi0) +abs(phi0));
    end
    iter_ct = 0; 
    swtch = 1;
    while abs(w) > erretm * eps
        %%% maximal iterations
        if iter_ct >= MAX_ITER
            if displaywarning
                warning('root %i does not converge after %i iteration, with residual %e\n', i, MAX_ITER, abs(w))
            end
            break;    
        end
        
        %%% update upper and lower bounds
        if w >= 0
            ub = min(ub, x0);
        else
            lb = max(lb, x0);
        end
        
        % swtch = 1 : fixed weight method
        % swtch = 0 : middle way method
        if iter_ct > 0
            if (prew * w > 0) && (abs(prew) > abs(w) / 10)
                swtch = - swtch + 1;
            end
        end
        
        %%% quadratic equation coefficients
        a = (delta(i) - x0 + delta(i+1) - x0) * w - (delta(i) - x0) * (delta(i+1) - x0) * dw;
        b = (delta(i) - x0) * (delta(i+1) - x0) * w;
        if swtch == 0   
            c =  - (delta(i) - x0) * dpsi(x0, delta, i) - (delta(i+1) - x0) * dphi(x0, delta, i) + w;       
        elseif swtch == 1 
            if f0(i) >= 0
                c = - (delta(i+1) - x0) * dw - v(i)^2 * (delta(i) - delta(i+1)) / (delta(i) - x0)^2 + w;  
            elseif f0(i) < 0
                c = - (delta(i) - x0) * dw - v(i+1)^2 * (delta(i+1) - delta(i)) / (delta(i+1) - x0)^2 + w;  
            end
        end
        
        %%% eta : root update
        if abs(c) == 0          % handle corner case : c == 0
            if abs(a) == 0
                if f0(i) >= 0
                    a = v(i)^2 + (delta(i+1) - x0)^2 * (dw - v(i)^2/(delta(i) - x0)^2);
                elseif f0(i) < 0
                    a = v(i+1)^2 + (delta(i) - x0)^2 * (dw - v(i+1)^2/(delta(i+1) - x0)^2);
                end
            end
            eta = b / a;
        else    
            if a <= 0 
                eta = (a - sqrt(abs(a^2 - 4*b*c))) / (2*c);
            else
                eta = 2*b / (a + sqrt(abs(a^2 - 4*b*c)));
            end
        end
        
        
        %%% f*eta should be negative, otherwise run a Newton step
        if w*eta >= 0
            eta = - w / dw;
        end
        
        
        %%% check if updated root lies in the [lb, ub]
        if ((x0 + eta) < lb) || ((x0 + eta) > ub) || ((x0 + eta) <= lb0) || ((x0 + eta) >= ub0)
            if w >= 0
                eta = (lb - x0) / 2;
            else
                eta = (ub - x0) / 2;
            end
        end
        
        
        %%% update root
        x0 = x0 + eta;
        
        
        %%% for next iteration
        iter_ct = iter_ct + 1;
        prew = w;
        psi0 = psi(x0, delta, i); 
        phi0 = phi(x0, delta, i);
        dpsi0 = dpsi(x0, delta, i); 
        dphi0 = dphi(x0, delta, i);
        w = rho + psi0 + phi0;
        dw = dpsi0 + dphi0;
        nflops = nflops + 16*n;
        if stop_criteria  == 0
            erretm = abs(dot([org0-1:-1:1]', v(1:org0-1).^2 ./ (delta(1:org0-1)-x0))) + ...
                                  dot([n-org0:-1:1]', v(n:-1:org0+1).^2 ./ (delta(n:-1:org0+1)-x0));
            erretm = 8 * (phi(x0, delta, org0) - psi(x0, delta, org0-1)) + erretm + 2 * rho...
                            + 3 * abs(v(org0)^2 / (delta(org0)-x0)) + abs(x0) * dw;
        elseif stop_criteria  == 1
            erretm = C * n * (rho + abs(psi0) +abs(phi0));
        end
    end
    if record
        residual(i) = abs(w);
        erretms(i) = abs(erretm);
        iter(i) = iter_ct;
    end

    % i^th root and eigenvector
    if org0 == i
        tau(i) = x0;
        x(i) = x0 + d(i);
    elseif org0 == i+1
        tau(i) = x0;
        x(i) = x0 + d(i+1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% find n^th root %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% set up
delta = d - d(n);
psi = @(x) dot(v(1:n-1).^2, 1./(delta(1:n-1) - x));
dpsi = @(x) dot(v(1:n-1).^2, 1./(delta(1:n-1) - x).^2);
phi = @(x) v(n)^2 / (delta(n) - x);
dphi = @(x) v(n)^2 / (delta(n) - x)^2;

%%% initial guess
x0 = 1 / (2 * rho);
lb0 = 0;
w = rho + psi(x0) + phi(x0);
c = w - v(n-1)^2 / (delta(n-1) - x0) - v(n)^2 / (delta(n) - x0);
if w <= 0
    temp = v(n-1)^2 / (1/rho - delta(n-1)) + v(n)^2 * rho;
    if c <= temp
        x0 = 1 / rho;
    else
        a = c * delta(n-1) + v(n-1)^2 + v(n)^2;
        b = - v(n)^2 * delta(n-1);
        if a < 0
            x0 = 2*b / (sqrt(abs(a^2+4*b*c)) - a);
        else
            x0 = (a + sqrt(abs(a^2+4*b*c))) / (2*c);
        end
    end
    lb = 1 / (2 * rho);
    ub = 1 / rho;
else
    a = c * delta(n-1) + v(n-1)^2 + v(n)^2;
    b = - v(n)^2 * delta(n-1);
    if a < 0
        x0 = 2*b / (sqrt(abs(a^2+4*b*c)) - a);
    else
        x0 = (a + sqrt(abs(a^2+4*b*c))) / (2*c);
    end
    lb = 0;
    ub = 1 / (2 * rho);
end

%%% rational interpolation
psi0 = psi(x0); 
phi0 = phi(x0);
dpsi0 = dpsi(x0); 
dphi0 = dphi(x0);
w = rho + psi0 + phi0;
dw = dpsi0 + dphi0;
nflops = nflops + 16*n;

if stop_criteria == 0
    erretm = abs(dot([n-1:-1:1]', v(1:n-1).^2 ./ (delta(1:n-1)-x0)));
    erretm = 8 * (- psi0 - phi0) + erretm - phi0 + rho + abs(x0) * dw;
elseif stop_criteria == 1
    erretm = C * n * (rho + abs(psi0) + abs(phi0));
end

iter_ct = 0;
while abs(w) > erretm * eps 
    %%% maximal iterations
    if iter_ct >= MAX_ITER
        if displaywarning
            warning('root %i does not converge after %i iteration, with residual %e\n', n, MAX_ITER, abs(w))
        end
        break
    end
    
    
    %%% update upper and lower bounds
    if w >= 0
        ub = min(ub, x0);   
    else
        lb = max(lb, x0);
    end

    
    %%% calculate new root
    a = (delta(n-1) - x0 + delta(n) - x0) * w - (delta(n-1) - x0 ) * (delta(n) - x0) * dw;
    b = (delta(n-1) - x0 ) * (delta(n) - x0) * w;
    c = - (delta(n-1) - x0) * dpsi(x0) - (delta(n) - x0) * dphi(x0) + w;
    c = abs(c);
    if c == 0
        eta = ub - x0;
    else
        if a >= 0 
            eta = (a + sqrt(abs(a^2 - 4*b*c))) / (2*c);
        else
            eta = 2*b / (a - sqrt(abs(a^2 - 4*b*c)));
        end
    end
    
    %%% f*eta should be negative, otherwise run a Newton step
    if w*eta >= 0
        eta = -w / dw;
    end

    
    %%% check if updated root lies in the [lb, ub]
    if ((x0 + eta) < lb) || ((x0 + eta) > ub) || ((x0 + eta) <= lb0)
        if w >= 0
            eta = (lb - x0) / 2;
        else
            eta = (ub - x0) / 2;
        end
    end
    
    
    %%% update root
    x0 = x0 + eta;

    %%% for next iteration
    psi0 = psi(x0); 
    phi0 = phi(x0);
    dpsi0 = dpsi(x0); 
    dphi0 = dphi(x0);
    w = rho + psi0 + phi0;
    dw = dpsi0 + dphi0;
    nflops = nflops + 16*n;
    
    if stop_criteria == 0
        erretm  = abs(dot([n-1:-1:1]', v(1:n-1).^2 ./ (delta(1:n-1)-x0)));
        erretm = 8 * (- psi0 - phi0) + erretm - phi0 + rho + abs(x0) * dw;
    elseif stop_criteria == 1
        erretm = C * n * (rho + abs(psi0) + abs(phi0));
    end
    
    iter_ct = iter_ct + 1;
end

if record
    residual(n) = abs(w);
    erretms(n) = abs(erretm);
    iter(n) = iter_ct;
end

%%% n^th root and eigenvector
tau = [tau; x0];
org = [org; n];
x = [x; x0 + d(n)];






