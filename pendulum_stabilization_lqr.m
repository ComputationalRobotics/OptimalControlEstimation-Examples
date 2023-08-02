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

Q = eye(2); 
R = 1;

[K,~,~] = dlqr(A,B,Q,R,zeros(2,1));

z0 = [0.1;0.1];
num_steps = 1000;

z_traj = zeros(2,num_steps+1);
z_traj(:,1) = z0;
u_traj = zeros(1,num_steps);

z = z0;
for i = 1:num_steps
    ui = -K*z;
%     ui = clip(-K*z,-5,5);
    zdoti = pendulum_f_z(z,ui,m,g,l,b);
    z = zdoti * dt + z;
    
    u_traj(:,i) = ui;
    z_traj(:,i+1) = z;
end

%% plot comparison
labelsize = 20;

t_traj = (1:(num_steps)) * dt;
figure;
tiledlayout(3,1)
nexttile
plot(t_traj,z_traj(1,1:end-1),'LineWidth',2)
ylabel('$z_1$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;

nexttile
plot(t_traj,z_traj(2,1:end-1),'LineWidth',2)
ylabel('$z_2$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;

nexttile
plot(t_traj,u_traj,'LineWidth',2)
ylabel('$u$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
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

function u_clip = clip(u,umin,umax)
u_clip = min(umax,max(u,umin));
end