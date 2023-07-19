clc; clear; close all;

sostoolspath = "../SOSTOOLS";
mosekpath    = "../mosek";
addpath(genpath(sostoolspath))
addpath(genpath(mosekpath))

m_1 = 1; m_2 = 1; l_1 = 1; l_2 = 1; l_c1 = l_1/2; l_c2 = l_2/2; g = 9.8;
I_1 = 1/3*m_1*l_1^2/3; I_2 = 1/3*m_2*l_2^2;
ns = 6; ny = 4;

gamma = 0.0;  % desired exponential rate
deg_V = 2; % degree of Lyapunov function V
deg_K = 2; % degree of observer gain
kappa_plus = 0; % choose whether to go above the minimum relaxation order

deg_Q = deg_V - 2;
deg_M = deg_Q + deg_K;

e = mpvar('e',[ns,1]);
x = mpvar('x',[ns,1]);
w = mpvar('w',[3,1]); % note that W is symmetric hence we only need 3 vars
C = [eye(ny), zeros(ny, ns - ny)];
Ce = C*e;
Cx = C*x;
delta_f = acrobot_f(x+e,w,m_2,l_1,l_c2) - acrobot_f(x,w,m_2,l_1,l_c2);
deg_delta_f = delta_f.maxdeg;

M_q = [I_1+I_2+m_2*l_1^2+2*m_2*l_1*l_c2*x(4), I_2+m_2*l_1*l_c2*x(4);
    I_2+m_2*l_1*l_c2*x(4), I_2];
W = [w(1), w(2); w(2), w(3)];
WM = W * M_q;

% equality constraints on x, e and H
h = [x(1)^2 + x(2)^2 - 1;
     x(3)^2 + x(4)^2 - 1;
     WM(:) - 1;];

% inequality constraints on x and e
g = [1;
     4 - e(1:4).^2;
     25 - e(5:6).^2;];

max_deg = max([deg_V-1+deg_delta_f, ...
    2+deg_M, ...
    h.maxdeg, ...
    g.maxdeg]);

% minimum relaxation order 
kappa   = ceil(max_deg / 2) + kappa_plus;

fprintf("deg_V: %d, deg_K: %d, kappa: %d.\n",deg_V,deg_K,kappa);

prog = sosprogram([e;x;w]);

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

eps = 0.001;

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

save('SOS-sols/acrobot_sol.mat', 'sol')

%% helper function
% State defined as x = [sin(theta_1), cos(theta_1), sin(theta_2),
% cos(theta_2), theta_1_dot, theta_2_dot]
function f = acrobot_f(x,t,m_2,l_1,l_c2)
    T = [t(1), t(2);
         t(2), t(3)];
    C = m_2*l_1*l_c2 * [-2*x(3)*x(6), -x(3)*x(6);
        x(3)*x(5), 0];
    f = [x(2)*x(5);
         -x(1)*x(5);
         x(4)*x(6);
         -x(3)*x(6);
    -T*C*[x(5); x(6)]];
end
