clc; clear; close all;

m = 1; l = 1; g = 9.8; b = 0.1; umax = 2;
initial_state = [0;0]; % start from the bottom
N = 50;
T0 = 6; % initial guess on T
x0 = zeros(2,N); % initial guess on (x1,...,x_N)
u0 = zeros(N,1); % initial guess on (u1,...,u_N)
x0(:,1) = initial_state;

% % bang-bang for initial guess
% for i = 1:N-1
%     u = -0.99*sign(x(2));
%     xdot = [x(2);-(b*x(2)+m*g*l*sin(x(1))+u)/m/l^2];
%     x = x + T0/(N-1)*xdot;
%     x0(:,i+1) = x;
%     u0(i+1,:) = u;
% end

v0 = [T0; x0(:); u0(:)];
lb = [0.1;
      -Inf*ones(2*N,1);
      -umax*ones(N,1)];
ub = [20;
      Inf*ones(2*N,1);
      umax*ones(N,1)];

obj = @(v) objective(v,N);
nonlincon = @(v) collocation(v,N,initial_state,m,l,g,b);

options = optimoptions('fmincon','Algorithm','interior-point',...
    'display','iter','MaxFunctionEvaluations',1e6,'MaxIterations',1e4);

[vopt,fopt,~,out] = fmincon(obj,v0,...
    [],[],[],[],... % no linear (in)equality constraints
    lb,ub,...
    nonlincon,options);

fprintf("Maximum constraint violation: %3.2f.\n",out.constrviolation);
fprintf("Objective: %3.2f.\n",fopt);

%% plot solution
uopt = vopt(2*N+2:end);
xopt = vopt(2:2*N+1);
xopt = reshape(xopt,2,N);
Topt = vopt(1);
t_grid = linspace(0,1,N)*Topt;
figure;
scatter(t_grid,uopt,100,'filled'); hold on;
plot(t_grid,uopt,'LineWidth',2);
xlabel('$t$','FontSize',24,'Interpreter','latex');
ylabel('$u(t)$','FontSize',24,'Interpreter','latex');
ax = gca; ax.FontSize = 20;
grid on;
figure;
scatter(t_grid,xopt,100,'filled'); hold on;
plot(t_grid,xopt,'LineWidth',2);
xlabel('$t$','FontSize',24,'Interpreter','latex');
ylabel('$x(t)$','FontSize',24,'Interpreter','latex');
grid on;
legend('$\theta(t)$','$\dot{\theta}(t)$','FontSize',20,'Interpreter','latex');


% % integrate from the initial state using the found controls
% tt = linspace(0,Topt,1000);
% [t,sol] = ode89(@(t,y) pendulum_ode(t,y,uopt,t_grid,m,l,g,b),tt,initial_state);
% figure;
% plot(tt,sol,'LineWidth',2);
% xlabel('$t$','FontSize',24,'Interpreter','latex');
% ylabel('$x(t)$','FontSize',24,'Interpreter','latex');
% ax = gca; ax.FontSize = 20;
% legend('$x_1(t)$','$x_2(t)$','FontSize',24,'Interpreter','latex');
% grid on;


%% helper functions
function f = objective(v,N)
T = v(1); h = T / (N-1);
x = v(2:2*N+1); x = reshape(x,2,N);
theta = x(1,:);
thetadot = x(2,:);
u = v(2*N+2:end);

tmp = [cos(theta);
       sin(theta);
       thetadot];

x_cost = sum((tmp - [-1;0;0]).^2,1);

u_cost = u.^2;

f = sum(x_cost(:) + u_cost(:)) * h;
end


function dx = pendulum(x,u,m,l,g,b)
    dx = [x(2);
          -1/(m*l^2) * (-u + b*x(2)+m*g*l*sin(x(1))) ];
end

function [c,ceq] = collocation(v,N,initial_state,m,l,g,b)
T = v(1);
h = T/(N-1);
x = reshape(v(2:2*N+1),2,N);
u = v(2*N+2:end);

c = [];
ceq = [];

for k=1:N-1
    uk = u(k);
    ukp1 = u(k+1);
    xk = x(:,k);
    xkp1 = x(:,k+1);
    fk = pendulum(xk,uk,m,l,g,b);
    fkp1 = pendulum(xkp1,ukp1,m,l,g,b);
    
    % collocation points
    xkc = 0.5*(xk+xkp1) + h/8 * (fk - fkp1);
    ukc = 0.5*(uk + ukp1);
    dxkc = -3/(2*h) * (xk-xkp1) - 0.25*(fk + fkp1);
    
    % collocation constraint
    ceq = [ceq;
           dxkc - pendulum(xkc,ukc,m,l,g,b)];
end

ceq = [ceq;
       x(:,1) - initial_state; % initial condition
       x(:,end) - [pi;0]]; % land at target point
end

function dx = pendulum_ode(t,states,u_grid,t_grid,m,l,g,b)
u_t = interp1(t_grid,u_grid,t); % piece-wise linear
x = states;
    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = -(b*x(2)+m*g*l*sin(x(1))+u_t)/m/l^2;
end