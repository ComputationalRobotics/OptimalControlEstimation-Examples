clc; clear; close all;

m = 1; g = 9.8; l = 1; b = 0.1;

Ac = [
    0, 1;
    g/l, -b / (m*l^2)
];

Bc = [0; 1/(m*l^2)];

dt = 0.01;

A = dt*Ac + eye(2);
B = dt*Bc;

Q = diag([1,1]); 
R = 1;
C = [1, 0];
sigma = 0.02;
M = diag([0,sigma^2]);
N = sigma^2;

[K,~,~] = dlqr(A,B,Q,R,zeros(2,1));
[L,SIG,~] = dlqr(A,C',M,N,zeros(2,1));
L = L';
% L = SIG*C'*inv(C*SIG*C'+N);

num_steps = 1000;
x = [2;-2];
xhat = [0;0];

xhat_traj = zeros(2,num_steps+1);
x_traj = zeros(2,num_steps+1);
x_traj(:,1) = x;
xhat_traj(:,1) = xhat;
u_traj = zeros(1,num_steps);

for k = 1:num_steps
    % design control
    uk = -K*xhat;
    u_traj(k) = uk;

    % update true state
    % x = A*x + B*uk + [0; sigma*randn];
    xdot = pendulum_f_z(x,uk,m,g,l,b);
    x = xdot * dt + x + [0; sigma*randn];
    x_traj(:,k+1) = x;

    % generate measurement
    y = C*x + sigma*randn;

    % update xhat
    xhat = A*xhat + B*uk + L*(y - C*(A*xhat + B*uk));
    xhat_traj(:,k+1) = xhat;
end

%% plot comparison
labelsize = 20;

t_traj = (1:(num_steps)) * dt;
figure;
tiledlayout(2,1)
nexttile
plot(t_traj,x_traj(1,1:end-1),'LineWidth',2)
hold on
plot(t_traj,xhat_traj(1,1:end-1),'red','LineWidth',2)
ylabel('$x_1,\hat{x}_1$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
legend('$x_1$','$\hat{x}_1$','Interpreter','latex','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;

nexttile
plot(t_traj,x_traj(2,1:end-1),'LineWidth',2)
hold on
plot(t_traj,xhat_traj(2,1:end-1),'red','LineWidth',2)
ylabel('$x_2,\hat{x}_2$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
legend('$x_2$','$\hat{x}_2$','Interpreter','latex','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;



function zdot = pendulum_f_z(z,u,m,g,l,b)
z1 = z(1);
z2 = z(2);

zdot = [
    z2;
    1/(m*l^2)*(u - b*z2 + m*g*l*sin(z1))
];
end