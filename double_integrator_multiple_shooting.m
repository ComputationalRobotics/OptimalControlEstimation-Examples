clc; clear; close all;

initial_state = [-10;-0];
N = 51; % number of knot points
% add initial guess
v0 = [1; % guess on terminal time
      rand(2*N,1); % guess on states
      randn(N,1)]; % guess on control

lb = [0.01;
    -Inf*ones(2*N,1);
    -1*ones(N,1)];
ub = [10;
    Inf*ones(2*N,1);
    1*ones(N,1)];

obj = @(v) v(1); % minimize time
nonlincon = @(v) double_integrator_nonlincon(v,initial_state);

options = optimoptions('fmincon','Algorithm','interior-point',...
    'display','iter','MaxFunctionEvaluations',1e6,'MaxIterations',1e4);

[vopt,fopt,~,out] = fmincon(obj,v0,...
    [],[],[],[],... % no linear (in)equality constraints
    lb,ub,...
    nonlincon,options);

fprintf("Maximum constraint violation: %3.2f.\n",out.constrviolation);
fprintf("Objective: %3.2f.\n",fopt);

%% plot solution
uopt = vopt(2+2*N:3*N+1);
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
[t,sol] = ode45(@(t,y) double_integrator(t,y,vopt),tt,initial_state);
figure;
plot(tt,sol,'LineWidth',2);
xlabel('$t$','FontSize',24,'Interpreter','latex');
ylabel('$x(t)$','FontSize',24,'Interpreter','latex');
ax = gca; ax.FontSize = 20;
legend('$x_1(t)$','$x_2(t)$','FontSize',24,'Interpreter','latex');
grid on;


%% helper functions
function dx = double_integrator(t,states,v)
% return xdot at the selected times t and states, using information from v
% assume the controls in v define piece-wise constant control signal

T = v(1); % final time
N = (length(v) - 1) / 3; % number of knot points

u_grid = v(2+2*N:3*N+1); % N controls
t_grid = linspace(0,1,N)*T;

u_t = interp1(t_grid,u_grid,t,'previous'); % piece-wise constant

A = [0 1; 0 0];
B = [0; 1];

dx = A * states + B * u_t;
end

function [c,ceq] = double_integrator_nonlincon(v,initial_state)
% enforce x_k+1 = RK45(x_k, u_k); integration done using ode45

T = v(1); % final time
N = (length(v) - 1) / 3; % number of knot points
t_grid = linspace(0,1,N)*T;
x1 = v(2:N+1); % position 
x2 = v(2+N:2*N+1); % velocity
% u = v(2+2*N:3*N+1); % controls

% no inequality constraints
c = []; 

% equality constraints
ceq = [];
for i = 1:(N-1)
    ti = t_grid(i);
    tip1 = t_grid(i+1);

    xi = [x1(i);x2(i)];
    xip1 = [x1(i+1);x2(i+1)];

    % integrate system dynamics starting from xi in [ti,tip1]
    tt = ti:(tip1-ti)/20:tip1; % fine-grained time discretization
    [~,sol_int] = ode45(@(t,y) double_integrator(t,y,v),tt,xi);
    xip1_int = sol_int(end,:);

    % enforce them to be the same
    ceq = [ceq;
           xip1_int(1) - xip1(1);
           xip1_int(2) - xip1(2)];
end

% add initial state constraint
ceq = [ceq;
       x1(1) - initial_state(1);
       x2(1) - initial_state(2)];

% add terminal state constraint: land at origin
ceq = [ceq;
       x1(end);
       x2(end)]; 
end