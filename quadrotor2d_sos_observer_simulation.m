clc; clear; close all;

% Load solution from SOS program
sol = load('SOS-sols/quadrotor2d_sol.mat').sol;

sostoolspath = "../SOSTOOLS";
addpath(genpath(sostoolspath))

% Constants
g = 9.8;

% Model parameters
m = 1; l = 1; r = l/2; I = m*l^2/12;
ns = 7; ny = 4; nu = 2;
obs_idx = [1,3,5,6];
C = [
    1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 1, 0;
];

num_steps = 19990;
dt = 0.001;


%% simulate real dynamics
x0 = [1;0;1;0;1/2;1/2*sqrt(3);0];
u_traj  = zeros(nu,num_steps);
x_traj = zeros(ns,num_steps+1);
x_traj(:,1) = x0;
x = x0;
for i = 1:num_steps
    u1   = 1 + rand;
    u2   = 1 + rand;
    % u    = [u1;u2];
    u    = zeros(2,1);
    u_traj(:,i) = u;
    y    = C*x;
    xdot = quadrotor2d_f(x) + quadrotor2d_psi(u,y,m,g,I,r);

    x    = xdot * dt + x;
    x(5:6) = normc(x(5:6));
    x_traj(:,i+1) = x;
end
y_traj = x_traj(obs_idx,:);
theta_traj = atan2(x_traj(5,:),x_traj(6,:));


%% simulate observer
% xhat0 = [0;0;0;0;0;1;0.5];
xhat0 = x0;
xhat_traj = zeros(ns,num_steps+1);
xhat_traj(:,1) = xhat0;
xhat = xhat0;
Ces = [zeros(4,1)];
for i = 1:num_steps
    yi = y_traj(:,i);
    ui = u_traj(:,i);
    yhati = C*xhat;
    Ce = yhati - yi;    
    Ces = [Ces, Ce];
    Qi = Q_func(Ce,sol,ns,obs_idx);
    Mi = M_func(Ce,yi,sol,ns,obs_idx);

    xhatdot = quadrotor2d_f(xhat) + quadrotor2d_psi(u,yi,m,g,I,r) + (Qi \ Mi)*Ce;
    xhat = xhat + xhatdot * dt;

    xhat(5:6) = normc(xhat(5:6));
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


%% helper functions
function Q = Q_func(Ce,sol,ns,obs_idx)
    e = mpvar('e',[ns,1]);
    Q_tilde = double(subs(sol.Q, e(obs_idx), Ce));
    Q = Q_tilde + sol.eps*eye(ns);
end

function M = M_func(Ce,y,sol,ns,obs_idx)
    e = mpvar('e',[ns,1]);
    x = mpvar('x',[ns,1]);
    M = double(subs(sol.M, [e(obs_idx); x(obs_idx)], [Ce;y]));
end

% State defined as x = [x; xdot; y; ydot; sin(theta); cos(theta); theta_dot]
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

function psi = quadrotor2d_psi(u,y,m,g,I,r)
    psi =   [0;
            -(u(1)+u(2))*y(3)/m;
            0;
            (u(1)+u(2))*y(4)/m - g;
            0;
            0;
            1/I*r*(u(1)-u(2));
    ];
end
