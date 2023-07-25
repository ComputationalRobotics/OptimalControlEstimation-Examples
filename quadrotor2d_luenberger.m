clear; clc; close all; format compact;

cvxpath = "./cvx";
addpath(genpath(cvxpath))

% Constants
g = 9.8;

% Model parameters
m = 1; l = 1; r = l/2; I = m*l^2/12;
ns = 6; ny = 3; nu = 2;

% Define system matrices
A = zeros(ns,ns); A(1,2) = 1; A(3,4) = 1; A(5,6) = 1;
% psi defined at the bottom
C = zeros(ny,ns); C(1,1) = 1; C(2,3) = 1; C(3,5) = 1;
obs_idx = [1,3,5];

Ob = obsv(A,C);
if rank(Ob) ~= length(A)
    disp("(A,C) not observable!")
else
    disp("(A,C) is observable!")
end

% Find optimal K that maximizes convergence rate
cvx_begin
    cvx_solver mosek
    variable H(ns,ny)
    variable P(ns,ns) symmetric
    minimize(0)
    subject to
        P == semidefinite(ns)
        A'*P - C'*H' + P*A - H*C == - eye(ns);
cvx_end
max_gam = 0.5 / max(eig(P));
K = P \ H;


%% simulate real dynamics
dt = 0.01;
T = 50;
num_steps = floor(T/dt);

x0 = [1;0;1;0;0;0];
u_traj  = zeros(nu,num_steps);
x_traj = zeros(ns,num_steps+1);
x_traj(:,1) = x0;
x = x0;
for i = 1:num_steps
    u1   = 1 + rand;
    u2   = 1 + rand;
    u    = [u1;u2];
    u    = zeros(2,1);
    u_traj(:,i) = u;
    y    = C*x;
    xdot = A*x + quadrotor2d_psi(u,y,m,g,I,r);

    x    = xdot * dt + x;
    x_traj(:,i+1) = x;
end
y_traj = x_traj(obs_idx,:);


%% simulate observer
xhat0 = [0;0;0;0;0;0];
xhat_traj = zeros(ns,num_steps+1);
xhat_traj(:,1) = xhat0;
xhat = xhat0;
Ces = [zeros(4,1)];
for i = 1:num_steps
    yi = y_traj(:,i);
    ui = u_traj(:,i);
    % yhati = C*xhat;
    % Ce = yhati - yi;

    xhatdot = A*xhat + quadrotor2d_psi(u,yi,m,g,I,r) + K*C*(x-xhat);
    xhat = xhat + xhatdot * dt;

    xhat_traj(:,i+1) = xhat;
end

%% plot comparison
labelsize = 20;
e_norm_traj = sqrt(sum((xhat_traj - x_traj).^2,1));

t_traj = (1:(num_steps+1)) * dt;
figure;
tiledlayout(4,1)
nexttile
s_comp = [x_traj(5,:);xhat_traj(1,:)];
plot(t_traj,s_comp','LineWidth',2)
ylabel('$\sin \theta$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

nexttile
c_comp = [x_traj(6,:);xhat_traj(2,:)];
plot(t_traj,c_comp','LineWidth',2)
ylabel('$\cos \theta$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

nexttile
thdot_comp = [x_traj(7,:);xhat_traj(3,:)];
plot(t_traj,thdot_comp','LineWidth',2)
ylabel('$\dot{\theta}$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

nexttile
plot(t_traj,e_norm_traj','LineWidth',2)
ylabel('$\Vert \hat{x} - x \Vert$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

%% Helper functions
% State defined as x = [x; xdot; y; ydot; theta; theta_dot]
function psi = quadrotor2d_psi(u,y,m,g,I,r)
    psi =   [0;
            -(u(1,1)+u(2,1))*sin(y(3))/m;
            0;
            (u(1)+u(2))*cos(y(3))/m - g;
            0;
            1/I*r*(u(1)-u(2));
    ];
end

