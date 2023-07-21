clc; clear; close all;

m = 1; g = 9.8; l = 1; b = 0.1;
cvxpath = "../cvx";
addpath(genpath(cvxpath))

L = 200;
k = [1;1];
num_steps = 599999;
dt = 0.00001;
A = [-k(1) 1;-k(2) 0];
%A = [1 0;0 1];
%% confirm lambda 
% lambda = 0.249;
cvx_begin sdp
    variable P(2,2) symmetric
    variable lam(1,1)
    maximize(lam)
    subject to
        % diag(P) == 1
        A'*P+P*A+4*lam*P <=0
        P >= 0
cvx_end
lambda_max(P)
lambda_min(P)
lambda_max(A'*P  + P*A + 4*lambda*P)
lambda_min(A'*P  + P*A + 4*lambda*P)
coff_1 = -(2*lambda*L-b*sqrt(2/(L^2-1))*L*sqrt(lambda_max(P)/lambda_min(P)));
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
    xhatdot = [xhat(2)-L*k(1)*Ce;-(b*xhat(2)+m*g*l^2*sin(yi)-ui)/m/l^2-L^2*k(2)*Ce];
    xhat = xhat + xhatdot * dt;

    xhat_traj(:,i+1) = xhat;
end


%% plot comparison

coff_2 = log(lambda_max(P)/lambda_min(P))/2+log((xhat0(1)-x0(1))^2+((xhat0(2)-x0(2))/L)^2)/2;
labelsize = 20;
e_norm_traj = sqrt(sum((xhat_traj(1,:) - x_traj(1,:)).^2+(xhat_traj(2,:) - x_traj(2,:)).^2/L^2,1));

t_traj = (1:(num_steps+1)) * dt;
figure;
% tiledlayout(3,1)
% nexttile
% s_comp = [x_traj(1,:);xhat_traj(1,:)];
% plot(t_traj,s_comp','LineWidth',2)
% ylabel('$\theta$','Interpreter','latex','FontSize',labelsize)
% xlabel('time','FontSize',labelsize)
% ax = gca;
% ax.FontSize = 16;
% 
% nexttile
% c_comp = [x_traj(2,:);xhat_traj(2,:)];
% plot(t_traj,c_comp','LineWidth',2)
% ylabel('$\dot{\theta}$','Interpreter','latex','FontSize',labelsize)
% xlabel('time','FontSize',labelsize)
% ax = gca;
% ax.FontSize = 16;


nexttile
bound = coff_2+coff_1*t_traj;
plot(t_traj,log(e_norm_traj'),'LineWidth',2)
hold on;
plot(t_traj,bound','LineWidth',2)
ylabel('$log (\Vert \hat{x} - x \Vert)$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;
