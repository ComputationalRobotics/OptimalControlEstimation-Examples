clc; clear; close all;

sostoolspath = "../SOSTOOLS";
mosekpath    = "../../../mosek";
addpath(genpath(sostoolspath))
addpath(genpath(mosekpath))

m = 1; g = 9.8; l = 1; b = 0.1;

gamma = 0.2; % desired exponential rate
deg_V = 2; % user-chosen degree of V
deg_K = 2; % user-chosen degree of K
kappa_plus = 0; % if >0 then go above the minimum relaxation order

% compute other degrees
deg_Q = deg_V - 2;
deg_M = deg_Q + deg_K;

e = mpvar('e',[3,1]);
x = mpvar('x',[3,1]);
C = [1, 0, 0; 0, 1, 0];
Ce = C*e;
Cx = C*x;
h = x(1)^2 + x(2)^2 - 1;
deg_h = h.maxdeg;
eps = 0.1;

delta_f = pendulum_f(x+e,m,b,l) - pendulum_f(x,m,b,l);
deg_delta_f = delta_f.maxdeg;

max_deg = max([deg_V - 1 + deg_delta_f,2+deg_M,deg_h]);
% minimum relaxation order 
kappa   = ceil(max_deg / 2) + kappa_plus;

fprintf("deg_V: %d, deg_K: %d, kappa: %d.\n",deg_V,deg_K,kappa);

prog = sosprogram([e;x]);

[prog,V] = sossosvar(prog,monomials(e,0:deg_V));
[prog,Q] = sospolymatrixvar(prog,monomials(Ce,0:deg_Q),[3,3],'symmetric');
[prog,M] = sospolymatrixvar(prog,monomials([Ce;Cx],0:deg_M),[3,2]);

[prog,sig0] = sossosvar(prog,monomials([e;x],0:kappa));
[prog,lam1] = sospolyvar(prog,monomials([e;x],0:(2*kappa-deg_h)));

eq = e'*(Q+eps*eye(3))*delta_f + e'*M*(Ce) + ...
    sig0 + lam1 * h + gamma * V;

prog = soseq(prog,eq);
prog = sosmatrixineq(prog,Q);
prog = soseq(prog,diff(V,e)-e'*(Q+eps*eye(3)));
prog = soseq(prog,subs(V,e,zeros(3,1)));

options.solver = 'mosek';
prog = sossolve(prog,options);

threshold = 1e-6;
sol.V = cleanpoly(sosgetsol(prog,V),threshold);
sol.Q = cleanpoly(sosgetsol(prog,Q),threshold);
sol.M = cleanpoly(sosgetsol(prog,M),threshold);

save('SOS-sols/pendulum_sol.mat', 'sol')

%% helper function
function f = pendulum_f(x,m,b,l)
const = b / (m * l^2);
f = [x(2)*x(3);
    -x(1)*x(3);
    -const*x(3)
    ];
end

