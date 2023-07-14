clc; clear; close all;

sostoolspath = "../SOSTOOLS";
mosekpath    = "../../../mosek";
addpath(genpath(sostoolspath))
addpath(genpath(mosekpath))

% Model parameters
m = 1; l = 1; r = l/2; I = m*l^2/12;
ns = 7; ny = 4; nu = 2;

gamma = 0.0000001; % desired exponential rate
deg_V = 2; % user-chosen degree of V
deg_K = 2; % user-chosen degree of K
kappa_plus = 0; % if >0 then go above the minimum relaxation order

% compute other degrees
deg_Q = deg_V - 2;
deg_M = deg_Q + deg_K;

e = mpvar('e',[ns,1]);
x = mpvar('x',[ns,1]);
C = [
    1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 1, 0;
];
Ce = C*e;
Cx = C*x;
h = x(5)^2 + x(6)^2 - 1;
deg_h = h.maxdeg;
eps = 0.1;

delta_f = quadrotor2d_f(x+e) - quadrotor2d_f(x);
deg_delta_f = delta_f.maxdeg;

max_deg = max([deg_V - 1 + deg_delta_f,2+deg_M,deg_h]);
% minimum relaxation order 
kappa   = ceil(max_deg / 2) + kappa_plus;

fprintf("deg_V: %d, deg_K: %d, kappa: %d.\n",deg_V,deg_K,kappa);

prog = sosprogram([e;x]);

[prog,V] = sossosvar(prog,monomials(e,0:deg_V));
[prog,Q] = sospolymatrixvar(prog,monomials(Ce,0:deg_Q),[ns,ns],'symmetric');
[prog,M] = sospolymatrixvar(prog,monomials([Ce;Cx],0:deg_M),[ns,ny]);

[prog,sig0] = sossosvar(prog,monomials([e;x],0:kappa));
[prog,lam1] = sospolyvar(prog,monomials([e;x],0:(2*kappa-deg_h)));

eq = e'*(Q+eps*eye(ns))*delta_f + e'*M*(Ce) + ...
    sig0 + lam1 * h + gamma * V;

prog = soseq(prog,eq);
prog = sosmatrixineq(prog,Q);
prog = soseq(prog,diff(V,e)-e'*(Q+eps*eye(ns)));
prog = soseq(prog,subs(V,e,zeros(ns,1)));

options.solver = 'mosek';
prog = sossolve(prog,options);

threshold = 1e-4;
sol.V = cleanpoly(sosgetsol(prog,V),threshold);
sol.Q = cleanpoly(sosgetsol(prog,Q),threshold);
sol.M = cleanpoly(sosgetsol(prog,M),threshold);
sol.eps = eps;

save('SOS-sols/quadrotor2d_sol.mat', 'sol')

%% helper function
% State defined as x = [x, x_dot, y, y_dot, sin(theta), cos(theta), theta_dot]
function f = quadrotor2d_f(x)
f = [x(2);
    0;
    x(4);
    0;
    x(6)*x(7);
    -x(5)*x(7);
    0;
];
end
