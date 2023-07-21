clc; clear; close all;

mp = 1; g = 9.8; l = 1; b = 0.1; mc = 1;

L = 1;
k = [1;1];
num_steps = 59999;
dt = 0.001;


%% simulate real dynamics
x0 = [0;1;1;0];
u_traj  = zeros(1,num_steps);
x_traj = zeros(4,num_steps+1);
x_traj(:,1) = x0;
x = x0;
for i = 1:num_steps
    ui   = -1;
    if(x(2)<0)
        ui = -ui;
    end
    u_traj(i) = ui;
    y    = x(1:2);
    xdot = [x(3);x(4);(ui+mp*sin(x(2))*(l*x(4)^2+g*cos(x(2))))/(mc+mp*sin(x(2))^2);-(ui*cos(x(2))+mp*l*x(4)^2*cos(x(2))*sin(x(2))+(mp+mc)*g*sin(x(2)))/l/(mc+mp*sin(x(2))^2)];
    x    = xdot * dt + x;
    x_traj(:,i+1) = x;
end
y_traj = x_traj(1:2,:);


%% simulate observer
xhat0 = [0;1;4;5];
xhat_traj = zeros(4,num_steps+1);
xhat_traj(:,1) = xhat0;
xhat = xhat0;
for i = 1:num_steps
    yi = y_traj(:,i);
    ui = u_traj(i);
    yhati = xhat(1:2);
    Ce = yhati - yi;
    % not use y as yhat(1)
    %xhatdot = [xhat(2)-L*k(1)*Ce;-(b*xhat(2)+m*g*l^2*sin(xhat(1))-ui)/m/l^2-L^2*k(2)*Ce];
    %not do any thing
    xhatdot = [xhat(3)-L*k(1)*Ce(1);xhat(4)-L*k(1)*Ce(2);(ui+mp*sin(yi(2))*(l*xhat(4)^2+g*cos(yi(2))))/(mc+mp*sin(yi(2))^2)-L^2*k(2)*Ce(1);-(ui*cos(yi(2))+mp*l*xhat(4)^2*cos(yi(2))*sin(yi(2))+(mp+mc)*g*sin(yi(2)))/l/(mc+mp*sin(yi(2))^2)-L^2*k(2)*Ce(2)];
    %xhatdot = [0;0;0;0];
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
plot(t_traj,x_traj(3,:)','LineWidth',2)
hold on
plot(t_traj,xhat_traj(3,:)','LineWidth',2)
ylabel('$x$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
legend('x','xhat')
ax = gca;
ax.FontSize = 16;

nexttile
c_comp = [x_traj(4,:);xhat_traj(4,:)];
plot(t_traj,c_comp','LineWidth',2)
ylabel('$\theta$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;


nexttile
plot(t_traj,log(e_norm_traj'),'LineWidth',2)
ylabel('$log (\Vert \hat{x} - x \Vert)$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

