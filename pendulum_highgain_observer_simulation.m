clc; clear; close all;

m = 1; g = 9.8; l = 1; b = 0.1;

L = 10;
k = [1;1];
num_steps = 5999;
dt = 0.01;


%% simulate real dynamics
x0 = [1;0];
u_traj  = zeros(1,num_steps);
x_traj = zeros(2,num_steps+1);
x_traj(:,1) = x0;
x = x0;
for i = 1:num_steps
    ui   = 0;
    u_traj(i) = ui;
    y    = x(1);
    xdot = [x(2);-(b*x(2)+m*g*l^2*sin(x(1))-ui)/m/l^2];
    x    = xdot * dt + x;
    x_traj(:,i+1) = x;
end
y_traj = x_traj(1,:);


%% simulate observer
xhat0 = [1;3];
xhat_traj = zeros(2,num_steps+1);
xhat_traj(:,1) = xhat0;
xhat = xhat0;
for i = 1:num_steps
    yi = y_traj(:,i);
    ui = u_traj(i);
    yhati = xhat(1);
    Ce = yhati - yi;
    % not use y as yhat(1)
    xhatdot = [xhat(2)-L*k(1)*Ce;-(b*xhat(2)+m*g*l^2*sin(xhat(1))-ui)/m/l^2-L^2*k(2)*Ce];
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
s_comp = [x_traj(1,:);xhat_traj(1,:)];
plot(t_traj,s_comp','LineWidth',2)
ylabel('$\theta$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

nexttile
c_comp = [x_traj(2,:);xhat_traj(2,:)];
plot(t_traj,c_comp','LineWidth',2)
ylabel('$\dot{\theta}$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;


nexttile
plot(t_traj,log(e_norm_traj'),'LineWidth',2)
ylabel('$log (\Vert \hat{x} - x \Vert)$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

