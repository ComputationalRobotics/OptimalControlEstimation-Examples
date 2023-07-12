clc; clear; close all;

sostoolspath = "../SOSTOOLS";
mosekpath    = "../../../mosek";
addpath(genpath(sostoolspath))
addpath(genpath(mosekpath))

m_c = 1; m_p = 1; l = 1; g = 9.8;
e_bound = 100;
ns = 6;
ny = 4;

deg_Q = 2;
deg_V = 2;
deg_M = 2;
kappa = 2;
gamma = 0.1;  % TODO: update

e = mpvar('e',[ns,1]);
x = mpvar('x',[ns,1]);
C = [1, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 0, 1];
Ce = C*e;   
Cx = C*x;

prog = sosprogram([e;x]);

[prog,V] = sossosvar(prog,monomials(e,0:deg_V));
[prog,Q] = sospolymatrixvar(prog,monomials(Ce,0:deg_Q),[ns,ns],'symmetric');
[prog,M] = sospolymatrixvar(prog,monomials([Ce;Cx],0:deg_M),[ns,ny]);

h1 = x(3)^2 + x(4)^2 - 1;
deg_h1 = h1.maxdeg;
h2 = x(6)*m_c/m_p + x(6)*x(3)^2 - 1;
deg_h2 = h2.maxdeg;

g1 = - (e'*e - e_bound^2);
deg_g1 = g1.maxdeg;

[prog,sig0] = sossosvar(prog,monomials([e;x],0:kappa));
[prog,sig1] = sossosvar(prog,monomials([e;x],0:floor((2*kappa-deg_g1)/2)));
[prog,lam1] = sospolyvar(prog,monomials([e;x],0:(2*kappa-deg_h1)));
[prog,lam2] = sospolyvar(prog,monomials([e;x],0:(2*kappa-deg_h2)));

eq = e'*Q*(cartpole_f(x+e,l) - cartpole_f(x,l)) + e'*M*(Ce) + ...
    sig0 + lam1 * h1 + lam2 * h2 + gamma * V;

eps = 10^-2;
prog = soseq(prog,eq);
prog = sosmatrixineq(prog,Q);
prog = soseq(prog,diff(V,e)-e'*(Q+eps*eye(ns)));
prog = soseq(prog,subs(V,e,zeros(ns,1)));
options.solver = 'mosek';
prog = sossolve(prog,options);

threshold = 1e-6;
sol.V = cleanpoly(sosgetsol(prog,V),threshold);
sol.Q = cleanpoly(sosgetsol(prog,Q),threshold);
sol.M = cleanpoly(sosgetsol(prog,M),threshold);


%% helper function
% State defined as x = [x, x_dot, sin(theta), cos(theta), theta_dot, a]
% where a = 1/(m_c/m_p+sin(theta)^2)
function f = cartpole_f(x,l)
f = [x(2);
    x(6) * x(3)*l*x(5)^2;
    x(4)*x(5);
    -x(3)*x(5);
    -x(6)*x(5)^2*x(3)*x(4);
    -2*x(3)*x(4)*x(5)*x(6)^2;
    ];
end
