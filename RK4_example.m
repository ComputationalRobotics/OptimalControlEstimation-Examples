% time step and num of step
clc; clear; close all
h = 0.01;
N = 10000;
% initial x_0
x_ini = [1,0,0];

%% for one turn
% u_0...u_{k-1}
u0 = 0;
u_traj_fix = randn(1,N-1);
u_traj = [u0,u_traj_fix];
x_traj = RK4(@(x,u)f_example(x,u),u_traj,h,N,x_ini);
figure;
plot(0:h:N*h,x_traj(:,1)','g',0:h:N*h,x_traj(:,2)','b');
legend('x_1','x_2');
%% plot x_N as a function of u_0
M = 100;
x_N_traj = zeros(M,3);
epi = 1;
u_traj_fix = randn(1,N-1);
for i = 1:M
    u0 = i*epi;
    u_traj = u0*ones(N,1);
    %x_traj = RK4(@(x,u)f_example(x,u),u_traj,h,N,x_ini);
    x_traj = RK4(@(x,u)f_example_2(x,u),u_traj,h,N,x_ini);
    x_N_traj(i,:) = x_traj(N+1,:);
end
figure;
plot((epi:epi:epi*M),x_N_traj(:,1)','g',...
    (epi:epi:epi*M),x_N_traj(:,2)','b',...
    (epi:epi:epi*M),x_N_traj(:,3)','r',...
    'LineWidth',2);
xlabel('$u$','FontSize',20,'Interpreter','latex')
ylabel('$x$','FontSize',20,'Interpreter','latex')
legend('$x_1$','$x_2$','$x_3$','FontSize',20,'Interpreter','latex');
ax = gca;
ax.FontSize = 16;
grid on
%% help function
function y = RK4(f,u,h,N,xini)
    x = xini';
    y = zeros(N+1,size(x,1));
    y(1,:) = x';
    for i = 1:N
        k1 = f(x,u(i));
        k2 = f(x+h*k1/2,u(i));
        k3 = f(x+h*k2/2,u(i));
        k4 = f(x+h*k3,u(i));
        x = x + (k1+2*k2+2*k3+k4)/6*h; 
        y(i+1,:) = x';
    end
end

function y = f_example(x,u)
    y = zeros(2,1);
    y(1) = (1-x(2)^2)*x(1)-x(2)+u;
    y(2) = x(1);
end

function y = f_example_2(x,u)
    y = zeros(3,1);
    y(1) = 10*(x(2)-x(1));
    y(2) = x(1)*(u-x(3))-x(2);
    y(3) = x(1)*x(2)-3*x(3);
end