clc; clear; close all;

m = 1; g = 9.8; l = 1; b = 0.1;
cvxpath = "../cvx";
addpath(genpath(cvxpath))

L = 200;
k = [1;1];
num_steps = 5999;
dt = 0.001;
A = [-k(1) 1;-k(2) 0];
J = diag([1;1;1]);
J_inv = inv(J);
rng(1);
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
eulerangle0 = [0;0;0];
q0 = [1;0;0;0];
tr0 = [0;0;0];
v0 = [1;0;1];
w0 = [0;0;0];
x0 = [tr0;v0;w0;q0];
rotorpos1 = [1;0;0];
rotorpos2 = [0;1;0];
u_traj  = zeros(6,num_steps+1);
x_traj = zeros(13,num_steps+1);
rotor_traj = zeros(6,num_steps+1);
%seperate x and R because R is a matrix
x_traj(:,1) = x0;
rotor_traj(:,1) = [rotorpos1;rotorpos2]; 
x = x0;
for i = 1:num_steps
    %force we have in the body frame
    f = [0;0;9];
    %torque we have in the body frame
    tor = [50;-50;0];
    ui   = [f;tor];
    u_traj(:,i) = ui;
    R = (eye(3) + 2*x(10)*wHat_func(x(11:13)) + 2*wHat_func(x(11:13))*wHat_func(x(11:13)));
    xdot = [x(4:6);f/m-[0;0;g];-J_inv*cross(x(7:9),J*x(7:9))+tor;1/2*[-x(11:13)';x(10)*eye(3)+wHat_func(x(11:13))]*x(7:9)];
    % xdot = [x(4:6);(eye(3)+2*x(10)*wHat_func(x(11:13))+2*wHat_func(x(11:13))*wHat_func(x(11:13)))*f/m-[0;0;g];-J_inv*cross(x(7:9),J*x(7:9))+tor;1/2*[-x(11:13)';x(10)*eye(3)+wHat_func(x(11:13))]*x(7:9)];
    x    = xdot * dt + x;
    x(10:13) = normalize(x(10:13),'norm',2);
    % Rdot = R*crossmatrix(x(7:9));
    % R    = Rdot * dt + R;
    % R = expmap(R);
    x_traj(:,i+1) = x;
    rotor_traj(:,i+1) = [R*rotorpos1;R*rotorpos2]; 
    if(mod(i,50)==0)
        plot3(x_traj(1,i+1),x_traj(2,i+1),x_traj(3,i+1),'o');
        hold on;
        plot3(x_traj(1,i+1)+rotor_traj(1,i+1),x_traj(2,i+1)+rotor_traj(2,i+1),x_traj(3,i+1)+rotor_traj(3,i+1),'o');
        hold on;
        plot3(x_traj(1,i+1)-rotor_traj(1,i+1),x_traj(2,i+1)-rotor_traj(2,i+1),x_traj(3,i+1)-rotor_traj(3,i+1),'o');
        hold on;
        plot3(x_traj(1,i+1)+rotor_traj(4,i+1),x_traj(2,i+1)+rotor_traj(5,i+1),x_traj(3,i+1)+rotor_traj(6,i+1),'o');
        hold on;
        plot3(x_traj(1,i+1)-rotor_traj(4,i+1),x_traj(2,i+1)-rotor_traj(5,i+1),x_traj(3,i+1)-rotor_traj(6,i+1),'o');
        axis([-10 10 -10 10 -5 5]);
        mov(i/50) = getframe;
        clf
    end
end
% movie(mov,1);
%R is not in here
y_traj = [x_traj(1:3,:);x_traj(10:13,:)];
%% draw simulation
% plot3(x_traj(1,:),x_traj(2,:),x_traj(3,:));
% hold on;
% plot3(x_traj(1,:)+rotor_traj(1,:),x_traj(2,:)+rotor_traj(2,:),x_traj(3,:)+rotor_traj(3,:));
% hold on;
% plot3(x_traj(1,:)-rotor_traj(1,:),x_traj(2,:)-rotor_traj(2,:),x_traj(3,:)-rotor_traj(3,:));
% hold on;
% plot3(x_traj(1,:)+rotor_traj(4,:),x_traj(2,:)+rotor_traj(5,:),x_traj(3,:)+rotor_traj(6,:));
% hold on;
% plot3(x_traj(1,:)-rotor_traj(4,:),x_traj(2,:)-rotor_traj(5,:),x_traj(3,:)-rotor_traj(6,:));
% hold on;
% axis equal
% axis tight
%% simulate observer
qhat0 = [1;0;0;0];
trhat0 = [0;0;0];
vhat0 = [1;0;1];
what0 = [0;0;0];
xhat0 = [tr0;v0;w0;q0];
xhat_traj = zeros(13,num_steps+1);
rotor_traj = zeros(6,num_steps+1);
%seperate x and R because R is a matrix
x_traj(:,1) = x0;
rotor_traj(:,1) = [rotorpos1;rotorpos2]; 
x = x0;
for i = 1:num_steps
    %force we have in the body frame
    f = [0;0;9];
    %torque we have in the body frame
    tor = [5;-5;0];
    ui   = [f;tor];
    u_traj(:,i) = ui;
    R = (eye(3)+2*x(10)*wHat_func(x(11:13))+2*wHat_func(x(11:13))*wHat_func(x(11:13)));
    xdot = [x(4:6);f/m-[0;0;g];-J_inv*cross(x(7:9),J*x(7:9))+tor;1/2*[-x(11:13)';x(10)*eye(3)+wHat_func(x(11:13))]*x(7:9)];
    % xdot = [x(4:6);(eye(3)+2*x(10)*wHat_func(x(11:13))+2*wHat_func(x(11:13))*wHat_func(x(11:13)))*f/m-[0;0;g];-J_inv*cross(x(7:9),J*x(7:9))+tor;1/2*[-x(11:13)';x(10)*eye(3)+wHat_func(x(11:13))]*x(7:9)];
    x    = xdot * dt + x;
    x(10:13) = normalize(x(10:13),'norm',2);
    % Rdot = R*crossmatrix(x(7:9));
    % R    = Rdot * dt + R;
    % R = expmap(R);
    x_traj(:,i+1) = x;
    rotor_traj(:,i+1) = [R*rotorpos1;R*rotorpos2]; 
    if(mod(i,50)==0)
        plot3(x_traj(1,i+1),x_traj(2,i+1),x_traj(3,i+1),'o');
        hold on;
        plot3(x_traj(1,i+1)+rotor_traj(1,i+1),x_traj(2,i+1)+rotor_traj(2,i+1),x_traj(3,i+1)+rotor_traj(3,i+1),'o');
        hold on;
        plot3(x_traj(1,i+1)-rotor_traj(1,i+1),x_traj(2,i+1)-rotor_traj(2,i+1),x_traj(3,i+1)-rotor_traj(3,i+1),'o');
        hold on;
        plot3(x_traj(1,i+1)+rotor_traj(4,i+1),x_traj(2,i+1)+rotor_traj(5,i+1),x_traj(3,i+1)+rotor_traj(6,i+1),'o');
        hold on;
        plot3(x_traj(1,i+1)-rotor_traj(4,i+1),x_traj(2,i+1)-rotor_traj(5,i+1),x_traj(3,i+1)-rotor_traj(6,i+1),'o');
        axis([-10 10 -10 10 -5 5]);
        mov(i/50) = getframe;
        clf
    end
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

%% functions
function wHat = wHat_func(x)
    wHat = [0 -x(3) x(2);x(3) 0 -x(1);-x(2) x(1) 0];
end