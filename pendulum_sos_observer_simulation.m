clc; clear; close all;

m = 1; g = 9.8; l = 1; b = 0.1;

num_steps = 5999;
dt = 0.01;

%% simulate real dynamics
x0 = [0;1;0];
u_traj  = zeros(1,num_steps);
x_traj = zeros(3,num_steps+1);
x_traj(:,1) = x0;
x = x0;
for i = 1:num_steps
    ui   = -1 + 2*rand;
    u_traj(i) = ui;
    y    = x(1:2);
    xdot = pendulum_f(x,m,b,l) + pendulum_psi(ui,y,m,g,l);

    x    = xdot * dt + x;
    x(1:2) = normc(x(1:2));
    x_traj(:,i+1) = x;
end
y_traj = x_traj(1:2,:);


%% simulate observer
xhat0 = [0;1;2];
xhat_traj = zeros(3,num_steps+1);
xhat_traj(:,1) = xhat0;
xhat = xhat0;
for i = 1:num_steps
    yi = y_traj(:,i);
    ui = u_traj(i);
    yhati = xhat(1:2);
    Ce = yhati - yi;
    Qi = Q_func(Ce);
    Mi = M_func(Ce,yi);
    
    xhatdot = pendulum_f(xhat,m,b,l) + pendulum_psi(ui,yi,m,g,l) + (Qi \ Mi)*Ce;

    xhat = xhat + xhatdot * dt;

    xhat(1:2) = normc(xhat(1:2));
    xhat_traj(:,i+1) = xhat;
end


%% plot comparison
labelsize = 20;
e_norm_traj = sqrt(sum((xhat_traj - x_traj).^2,1));

t_traj = (1:(num_steps+1)) * dt;
figure;
tiledlayout(4,1)
nexttile
s_comp = [x_traj(1,:);xhat_traj(1,:)];
plot(t_traj,s_comp','LineWidth',2)
ylabel('$\sin \theta$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

nexttile
c_comp = [x_traj(2,:);xhat_traj(2,:)];
plot(t_traj,c_comp','LineWidth',2)
ylabel('$\cos \theta$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

nexttile
thdot_comp = [x_traj(3,:);xhat_traj(3,:)];
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
function Q = Q_func(Ce)
e1 = Ce(1);
e2 = Ce(2);
Q = [
    0.4603*e2^2 + 1.1909, -0.4603*e1*e2, 0;
    -0.4603*e1*e2, 0.4603*e1^2 + 1.1909, 0;
    0, 0, 1.8863
];

end


function M = M_func(Ce,y)
e1 = Ce(1);
e2 = Ce(2);
M = [
    -2.0878*e1^2-0.8667*e2^2-0.4588-0.4885, -0.8667*e1*e2;
    -0.8667*e1*e2, -0.8667*e1^2-2.0878*e2^2-0.4588-0.4885;
    -1.1909*y(2),1.1909*y(1)
    ];
end

function f = pendulum_f(x,m,b,l)
const = b / (m * l^2);
f = [x(2)*x(3);
    -x(1)*x(3);
    -const*x(3)
    ];
end

function psi = pendulum_psi(u,y,m,g,l)
psi = [
    0;
    0;
    (u - m*g*l*y(1))/(m*l^2)
];
end