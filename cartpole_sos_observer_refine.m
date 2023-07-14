clc; clear; close all;

sostoolspath = "../SOSTOOLS";
mosekpath    = "../../../mosek";
addpath(genpath(sostoolspath))
addpath(genpath(mosekpath))

% state x = [z,s,c,zdot,thetadot]
ns = 5;
ny = 3;

gamma = 0.0;  % desired exponential rate
deg_V = 2;
deg_K = 2; % degree of observer gain
deg_Q = deg_V - 2;
deg_M = deg_Q + deg_K;
kappa_plus = 0; % choose whether to go above the minimum relaxation order

e = mpvar('e',[ns,1]);
x = mpvar('x',[ns,1]);
Ce = e(1:3);
Cx = x(1:3);

% equality constraints on x and e
h = [x(2)^2 + x(3)^2 - 1];
g = [monomials([e;x],0)];

eps = 0.1;
prog = sosprogram([e;x]);
[prog,V] = sossosvar(prog,monomials(e,0:deg_V));
[prog,Q] = sospolymatrixvar(prog,monomials(Ce,0:deg_Q),[ns,ns],'symmetric');
Qtld = Q + eps*eye(ns);
[prog,M] = sospolymatrixvar(prog,monomials([Ce;Cx],0:deg_M),[ns,ny]);

scaled_Vdot = e'*Qtld*cartpole_delta_f(x,e) + e'*M*Ce;

max_deg = max([deg_V,...
    full(max(sum(scaled_Vdot.degmat,2))),...
    h.maxdeg,g.maxdeg]);
kappa = ceil(max_deg / 2) + kappa_plus;

fprintf("deg_K: %d, kappa: %d.\n",deg_K,kappa);

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

eq = scaled_Vdot + gamma * V + sigs'*g + lams'*h;

prog = soseq(prog,eq);
prog = sosmatrixineq(prog,Q);
prog = soseq(prog,diff(V,e)-e'*Qtld);
prog = soseq(prog,subs(V,e,zeros(ns,1)));
options.solver = 'mosek';
prog = sossolve(prog,options);

threshold = 1e-6;
sol.Q = cleanpoly(sosgetsol(prog,Q),threshold);
sol.eps = eps;

save('SOS-sols/cartpole_tv_sol.mat', 'sol')

%% helper function
% state x = [z,s,c,zdot,thetadot]
function f = cartpole_delta_f(x,e)
z = x(1);
s = x(2);
c = x(3);
zdot = x(4); zdothat = zdot + e(4);
thetadot = x(5); thetadothat = thetadot + e(5);

f = [
    (1+s^2)*(zdothat - zdot);
    (1+s^2)*c*(thetadothat - thetadot);
    -(1+s^2)*s*(thetadothat - thetadot);
    s*(thetadothat^2 - thetadot^2);
    -c*s*(thetadothat^2 - thetadot^2)
];
end
