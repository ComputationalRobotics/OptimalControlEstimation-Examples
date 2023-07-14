clc; clear; close all;

sostoolspath = "../SOSTOOLS";
mosekpath    = "../../../mosek";
addpath(genpath(sostoolspath))
addpath(genpath(mosekpath))

% state x = [z,s,c,zdot,thetadot]
m_c = 1; m_p = 1; l = 1; g = 9.8;
ns = 5; n1 = 3; n2 = 2;
ny = n1;

% Choose whether to fix Q
FIX_Q = true;
RHO = 0.00001; % V(e) \leq RHO when FIX_Q is true
CONST = (3 - sqrt(5))/2;

gamma = 0.0;  % desired exponential rate
deg_K = 2; % degree of observer gain
kappa_plus = 0; % choose whether to go above the minimum relaxation order
deg_M = deg_K;

e = mpvar('e',[ns,1]);
x = mpvar('x',[ns,1]);
x1 = x(1:3);
x2 = x(4:ns);
e1 = e(1:3);
e2 = e(4:ns);

% equality constraints on x and e
h = [x(2)^2 + x(3)^2 - 1];
g = [monomials([e;x],0)];

eps = 0.0;
prog = sosprogram([e;x]);

if FIX_Q
    Q = CONST * eye(n1);
else
    [prog,Q] = sospolymatrixvar(prog,monomials([e;x],0),[n1,n1],'symmetric');
end
Qtld = Q + eps*eye(n1);

[prog,M1] = sospolymatrixvar(prog,monomials([e1;x1],0:deg_M),[n1,ny]);
[prog,K2] = sospolymatrixvar(prog,monomials([e1;x1],0:deg_K),[n2,ny]);

V = e1'*Qtld*e1 + e2'*cartpole_M(x,m_p,m_c,l)*e2;
if FIX_Q
    g = [g; RHO - V; RHO/CONST - e'*e];
end

v_delta = cartpole_v(x+e) - cartpole_v(x);
c_delta = (cartpole_Cor(x,m_p,l) - cartpole_Cor(x+e,m_p,l))*(x2+e2);
Vdot = 2*e1'*Qtld*v_delta + 2*e1'*M1*e1 + 2*e2'*c_delta + 2*e2'*K2*e1;

max_deg = max([ ...
    full(max(sum(Vdot.degmat,2))), ...
    full(max(sum(V.degmat,2))), ...
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

eq = Vdot + gamma * V + sigs'*g + lams'*h;

prog = soseq(prog,eq);
if ~FIX_Q
    prog = sosmatrixineq(prog,Q);
end
options.solver = 'mosek';
prog = sossolve(prog,options);

threshold = 1e-6;
if ~FIX_Q
    sol.Q = cleanpoly(sosgetsol(prog,Q),threshold);
end
sol.M1 = cleanpoly(sosgetsol(prog,M1),threshold);
sol.eps = eps;

save('SOS-sols/cartpole_tv_sol.mat', 'sol')

%% helper function
% state x = [z,s,c,zdot,thetadot]
function v = cartpole_v(x)
v = [
    x(4);
    x(3)*x(5);
    -x(2)*x(5)];
end

function Cor = cartpole_Cor(x,m_p,l)
Cor = [0, -m_p*l*x(5)*x(2);
       0, 0];
end

function M = cartpole_M(x,m_p,m_c,l)
M = [m_c+m_p, m_p*l*x(3);
     m_p*l*x(3), m_p*l^2];
end

