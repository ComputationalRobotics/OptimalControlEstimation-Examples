clc; clear; close all;

initial_state = [-10;0];
N = 51;

T0 = 1; % initial guess on T
x0 = randn(2,N); % initial guess on (x1,...,x_N)
u0 = randn(N,1); % initial guess on (u1,...,u_N)
v0 = [T0; x0(:); u0(:)];

lb = [0.1;
      -Inf*ones(2*N,1);
      -1*ones(N,1)];
ub = [10;
      Inf*ones(2*N,1);
      1*ones(N,1)];

obj = @(v) v(1); % minimize time
nonlincon = @(v) collocation(v,N,initial_state);

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
Topt = vopt(1);
t_grid = linspace(0,1,N)*Topt;
figure;
scatter(t_grid,uopt,100,'filled'); hold on;
plot(t_grid,uopt,'LineWidth',2);
xlabel('$t$','FontSize',24,'Interpreter','latex');
ylabel('$u(t)$','FontSize',24,'Interpreter','latex');
ax = gca; ax.FontSize = 20;
grid on;

% integrate from the initial state using the found controls
tt = linspace(0,Topt,1000);
[t,sol] = ode45(@(t,y) double_integrator_ode(t,y,uopt,t_grid),tt,initial_state);
figure;
plot(tt,sol,'LineWidth',2);
xlabel('$t$','FontSize',24,'Interpreter','latex');
ylabel('$x(t)$','FontSize',24,'Interpreter','latex');
ax = gca; ax.FontSize = 20;
legend('$x_1(t)$','$x_2(t)$','FontSize',24,'Interpreter','latex');
grid on;


%% helper functions
function dx = double_integrator(x,u)
A = [0 1; 0 0];
B = [0; 1];
dx = A * x + B * u;
end

function [c,ceq] = collocation(v,N,initial_state)
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
    fk = double_integrator(xk,uk);
    fkp1 = double_integrator(xkp1,ukp1);
    
    % collocation points
    xkc = 0.5*(xk+xkp1) + h/8 * (fk - fkp1);
    ukc = 0.5*(uk + ukp1);
    dxkc = -3/(2*h) * (xk-xkp1) - 0.25*(fk + fkp1);
    
    % collocation constraint
    ceq = [ceq;
           dxkc - double_integrator(xkc,ukc)];
end

ceq = [ceq;
       x(:,1) - initial_state; % initial condition
       x(:,end)]; % land at zero
end

function dx = double_integrator_ode(t,states,u_grid,t_grid)
u_t = interp1(t_grid,u_grid,t); % piece-wise linear
A = [0 1; 0 0];
B = [0; 1];
dx = A * states + B * u_t;
end