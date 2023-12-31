clc; clear; close all;
restoredefaultpath;

%% parameters of the pendulum
m = 1; g = 9.8; l = 1; b = 0.1;
umax = 2;
A = [0, 1; g/l, -b/(m*l^2)];
B = [0; 1/(m*l^2)];
Q = diag([1,1]); 
R = 1;

%% TODO: Compute K and S from LQR
[K,S] = ;

%% Region to draw samples
x1_lb = -0.2*pi;
x1_ub = 0.2*pi;
x2_lb = -0.2*pi;
x2_ub = 0.2*pi;

%% Simulate pendulum with clipped control
pendulum_sat = @(t,x) pendulum_f(x,clip(-K*x,-umax,umax),m,g,l,b);
T = 100;
N = 1000;
stable = [];
unstable = [];
for i = 1:N
    x1 = x1_lb + rand*(x1_ub - x1_lb);
    x2 = x2_lb + rand*(x2_ub - x2_lb);
    x0 = [x1;x2];
    [t,y] = ode89(pendulum_sat,[0,T],x0);
    if norm(y(end,:)) < 1e-3
        stable = [stable,x0];
    else
        unstable = [unstable,x0];
    end
end
% plot points
figure;
scatter(stable(1,:),stable(2,:),80,'filled','Marker','o');
hold on
scatter(unstable(1,:),unstable(2,:),80,'filled','Marker','square');
xlabel('$x_1$','Interpreter','latex','FontSize',24);
ylabel('$x_2$','Interpreter','latex','FontSize',24);
ax = gca; ax.FontSize = 20;
hold on

%% certify ROA using SOS programming
sostoolspath = '../SOSTOOLS';
mosekpath = '../../../mosek';
addpath(genpath(sostoolspath));
addpath(genpath(mosekpath));

x = mpvar('x',[2,1]);
J = x'*S*x;
eps = 0.01;
rho = 1;
kappa = 4; % relaxation order

prog = sosprogram(x);
% case 1
Jdot1 = - jacobian(J,x)*pendulum_poly(x,umax,m,g,l,b) - eps*(x'*x);
set1 = [rho - J; -K*x - umax];
[prog, rhs1] = SOSonSet(prog,x,kappa,set1);
prog = soseq(prog, Jdot1 - rhs1);

% TODO: case 2
Jdot2 = ;
set2 = ;
[prog, rhs2] = SOSonSet(prog,x,kappa,set2);
prog = soseq(prog, Jdot2 - rhs2);

% TODO: case 3
Jdot3 = ;
set3 = ;
[prog, rhs3] = SOSonSet(prog,x,kappa,set3);
prog = soseq(prog, Jdot3 - rhs3);

options.solver = 'mosek';
prog = sossolve(prog, options);

e = ellipse(S/rho,zeros(2,1));
plot(e(1,:),e(2,:),'green','LineWidth',4);
legend('Stabilized','Non-Stabilized','Certified','FontSize',22);


%% helper function
function v = clip(u,a,b)
%% TODO: implement the clip function that clips u between [a,b]
end

function zdot = pendulum_f(z,u,m,g,l,b)
% pendulum true dynamics
z1 = z(1);
z2 = z(2);

zdot = [
    z2;
    1/(m*l^2)*(u - b*z2 + m*g*l*sin(z1))
];
end

function zdot = pendulum_poly(z,u,m,g,l,b)
%% TODO: implement pendulum approximate polynomial dynamics
end

function [prog,rhs] = SOSonSet(prog,var,kappa,set_ineq)
% a polynomial p is SOS on a set K defined by inequality constraints
% relaxation order kappa
v0 = monomials(var,0:kappa);
[prog, sig0] = sossosvar(prog,v0);
rhs = sig0;
for j = 1:length(set_ineq)
    gj = set_ineq(j);
    degj = floor((2*kappa - gj.maxdeg)/2);
    vj = monomials(var,0:degj);
    [prog, sigj] = sossosvar(prog,vj);
    rhs = rhs + sigj*gj;
end
end

function e = ellipse(A,center)
[U,~] = eig(A);
theta = linspace(0,2*pi,1000);
elipoid = zeros(2,1000);
for pt = 1:1000
    elipoid(:,pt) = U(:,1)*sin(theta(pt)) + U(:,2)*cos(theta(pt));
    elipoid(:,pt) = elipoid(:,pt)/sqrt(elipoid(:,pt)'*A*elipoid(:,pt));
end
e = elipoid + center;
end
