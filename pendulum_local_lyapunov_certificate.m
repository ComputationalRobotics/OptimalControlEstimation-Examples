clc; clear; close all;

sostoolspath = "../SOSTOOLS";
mosekpath    = "../../../mosek";
addpath(genpath(sostoolspath))
addpath(genpath(mosekpath))

m = 1; g = 9.8; l = 1; b = 0.1;

deg_V = 2; % degree of Lyapunov function

x = mpvar('x',[3,1]);
x0 = [0;1;0]; % equilibrium point
prog = sosprogram(x);
[prog,V] = sospolyvar(prog,monomials(x,0:deg_V));
eps1 = 0.01;
eps2 = 0.01;
fx = pendulum_f(x,m,b,l,g);
Vdot = diff(V,x)*fx;
deg_Vdot = max(full(sum(Vdot.degmat,2)));

kappa = ceil( max([deg_V/2,deg_Vdot/2]) );

h = [x(1)^2 + x(2)^2 - 1]; % equality constraint
g = [monomials(x,0);
    (pi/2)^2 - x(3)^2; % velocity bounded
    x(2) - 3/4]; % cos theta \geq lower_bound

%% enforce V \geq eps1 * (x-x0)'*(x-x0)
lams = [];
for i = 1:length(h)
    deg_hi = h(i).maxdeg;
    [prog,lami] = sospolyvar(prog,monomials(x,0:(2*kappa-deg_hi)));
    lams = [lams;lami];
end
sigs = [];
for i = 1:length(g)
    deg_gi = g(i).maxdeg;
    [prog,sigi] = sossosvar(prog,monomials(x,0:floor((2*kappa-deg_gi)/2)));
    sigs = [sigs;sigi];
end

eq = V - (eps1 * (x-x0)'*(x-x0)) - sigs'*g - lams'*h;
prog = soseq(prog,eq);

%% enforce Vdot \leq -eps1 * (x-x0)'*(x-x0)
lams_dot = [];
for i = 1:length(h)
    deg_hi = h(i).maxdeg;
    [prog,lami_dot] = sospolyvar(prog,monomials(x,0:(2*kappa-deg_hi)));
    lams_dot = [lams_dot;lami_dot];
end
sigs_dot = [];
for i = 1:length(g)
    deg_gi = g(i).maxdeg;
    [prog,sigi_dot] = sossosvar(prog,monomials(x,0:floor((2*kappa-deg_gi)/2)));
    sigs_dot = [sigs_dot;sigi_dot];
end
eq_dot = (- eps2 * (x-x0)'*(x-x0)) - Vdot - sigs_dot'*g - lams_dot'*h;
prog = soseq(prog,eq_dot);

%% enforce V(x0) = 0, Vdot(x0) = 0
prog = soseq(prog,subs(V,x,x0));
prog = soseq(prog,subs(Vdot,x,x0));

options.solver = 'mosek';
prog = sossolve(prog,options);

sol.V = cleanpoly(sosgetsol(prog,V),1e-8);
sol.Vdot = [diff(sol.V,x(1)), diff(sol.V,x(2)), diff(sol.V,x(3))] * fx;

%% plot V and Vdot
thetadot_range = linspace(-pi/2,pi/2,100);
theta_range = linspace(-acos(3/4),acos(3/4),100);
[theta, thetadot] = meshgrid(theta_range, thetadot_range);
V_val = zeros(100,100);
Vdot_val = zeros(100,100);

for i = 1:100
    for j = 1:100
        thetai = theta(i,j);
        thetadoti = thetadot(i,j);
        xi = pendulum_xi([thetai;thetadoti]);
        Vi = subs(sol.V,x,xi);
        Vdoti = subs(sol.Vdot,x,xi);
        
        V_val(i,j) = Vi;
        Vdot_val(i,j) = Vdoti;
    end
end

%% generate plots
figure;
surf(theta,thetadot,V_val)
xlabel('$\theta$','FontSize',20,'Interpreter','latex')
ylabel('$\dot{\theta}$','FontSize',20,'Interpreter','latex')
zlabel('$V(x)$','FontSize',20,'Interpreter','latex')

figure;
surf(theta,thetadot,Vdot_val)
xlabel('$\theta$','FontSize',20,'Interpreter','latex')
ylabel('$\dot{\theta}$','FontSize',20,'Interpreter','latex')
zlabel('$\dot{V}(x)$','FontSize',20,'Interpreter','latex')



%% pendulum dynamics (uncontrolled)
function f = pendulum_f(x,m,b,l,g)
s = x(1);
c = x(2);
thetadot = x(3);
f = [c*thetadot;
    -s*thetadot;
    -1/(m*l^2)*(b*thetadot + m*g*l*s)];
end

function x = pendulum_xi(x_org)
x = [sin(x_org(1)); cos(x_org(1)); x_org(2)];
end