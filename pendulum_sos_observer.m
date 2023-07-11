clc; clear; close all;

sostoolspath = "../SOSTOOLS";
mosekpath    = "../../../mosek";
addpath(genpath(sostoolspath))
addpath(genpath(mosekpath))

m = 1; g = 9.8; l = 1; b = 0.1;

deg_Q = 2;
deg_V = 2;
deg_M = 2;
kappa = 2;
gamma = 0.1;

e = mpvar('e',[3,1]);
x = mpvar('x',[3,1]);
C = [1, 0, 0; 0, 1, 0];
Ce = C*e;
Cx = C*x;

prog = sosprogram([e;x]);

[prog,V] = sossosvar(prog,monomials(e,0:deg_V));
[prog,Q] = sospolymatrixvar(prog,monomials(Ce,0:deg_Q),[3,3],'symmetric');
[prog,M] = sospolymatrixvar(prog,monomials([Ce;Cx],0:deg_M),[3,2]);

h = x(1)^2 + x(2)^2 - 1;
deg_h = h.maxdeg;
[prog,sig0] = sossosvar(prog,monomials([e;x],kappa));
[prog,lam1] = sospolyvar(prog,monomials([e;x],0:(2*kappa-deg_h)));

eq = e'*Q*(pendulum_f(x+e,m,b,l) - pendulum_f(x,m,b,l)) + e'*M*(Ce) + ...
    sig0 + lam1 * h + gamma * V;

prog = soseq(prog,eq);
prog = sosmatrixineq(prog,Q);
prog = soseq(prog,diff(V,e)-e'*Q);
prog = soseq(prog,subs(V,e,zeros(3,1)));
% prog = sossetobj(prog,0);
options.solver = 'mosek';
prog = sossolve(prog,options);

sol.V = sosgetsol(prog,V);
sol.Q = sosgetsol(prog,Q);
sol.M = sosgetsol(prog,M);



%% helper function
function f = pendulum_f(x,m,b,l)
const = b / (m * l^2);
f = [x(2)*x(3);
    -x(1)*x(3);
    -const*x(3)
    ];
end

