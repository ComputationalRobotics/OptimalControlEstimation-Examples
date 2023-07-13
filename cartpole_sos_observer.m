clc; clear; close all;

sostoolspath = "../SOSTOOLS";
mosekpath    = "../../../mosek";
addpath(genpath(sostoolspath))
addpath(genpath(mosekpath))

m_c = 1; m_p = 1; l = 1; g = 9.8;
ns = 6;
ny = 4;

gamma = 0.1;  % desired exponential rate
deg_V = 2; % degree of Lyapunov function V
deg_K = 2; % degree of observer gain
kappa_plus = 0; % choose whether to go above the minimum relaxation order
e_bound = 100; % error bound

deg_Q = deg_V - 2;
deg_M = deg_Q + deg_K;

e = mpvar('e',[ns,1]);
x = mpvar('x',[ns,1]);
C = [1, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 0, 1];
Ce = C*e;   
Cx = C*x;
delta_f = cartpole_f(x+e,l) - cartpole_f(x,l);
deg_delta_f = delta_f.maxdeg;

% equality constraints on x and e
h = [x(3)^2 + x(4)^2 - 1;
    x(6)*m_c/m_p + x(6)*x(3)^2 - 1];

% inequality constraints on x and e
% note that I add g0 = 1
g = [1;
     -e'*e+e_bound^2];

max_deg = max([deg_V-1+deg_delta_f, ...
    2+deg_M, ...
    h.maxdeg, ...
    g.maxdeg]);

% minimum relaxation order 
kappa   = ceil(max_deg / 2) + kappa_plus;

fprintf("deg_V: %d, deg_K: %d, kappa: %d.\n",deg_V,deg_K,kappa);

prog = sosprogram([e;x]);

[prog,V] = sossosvar(prog,monomials(e,0:deg_V));
[prog,Q] = sospolymatrixvar(prog,monomials(Ce,0:deg_Q),[ns,ns],'symmetric');
[prog,M] = sospolymatrixvar(prog,monomials([Ce;Cx],0:deg_M),[ns,ny]);


lams = [];
for i = 1:length(h)
    deg_hi = h(i).maxdeg;
    [prog,lami] = sospolyvar(prog,monomials([e;x],0:(2*kappa-deg_hi)));
    lams = [lams;lami];
end
sigs = [];
for i = 1:length(g)
    deg_gi = g(i).maxdeg;
    [prog,sigi] = sossosvar(prog,monomials([e;x],0:floor((2*kappa-deg_gi)/2)));
    sigs = [sigs;sigi];
end

eps = 0.1;

eq = e'*(Q+eps*eye(ns))*delta_f + e'*M*(Ce) + gamma * V + ...
    sigs'*g + lams'*h;

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

save('SOS-sols/cartpole_sol.mat', 'sol')

%% helper function
% State defined as x = [x, x_dot, sin(theta), cos(theta), theta_dot, a]
% where a = 1/(m_c/m_p+sin(theta)^2)
function f = cartpole_f(x,l)
f = [x(2);
    x(6)*x(3)*l*x(5)^2;
    x(4)*x(5);
    -x(3)*x(5);
    -x(6)*x(5)^2*x(3)*x(4);
    -2*x(3)*x(4)*x(5)*x(6)^2;
    ];
end
