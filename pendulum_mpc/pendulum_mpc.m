clc; clear; close all;

% pendulum parameters
m = 1; l = 1; g = 9.8; b = 0.1; umax = g/2;

% system noise
noise = 1.0;

initial_state = [0;0];
Ts = 0.1; % solve trajectory optimization every Ts time
N  = 51; % fixed horizon for the trajectory optimization problem
h = 0.1; % time interval for trajectory optimization
num_steps = 100;

xi = initial_state;

state_trajectory = [];
time_trajectory = [];
control_trajory = [];
uopt = randn(N,1);
xopt = randn(2,N);
for i = 1:num_steps
    [uopt,xopt] = PendulumTrajOpt(N,h,xi,m,l,g,b,umax,uopt,xopt);
    
    % only choose the first control
    ui = uopt(1);

    % apply the first control to the plant for Ts duration
    [t,sol] = ode89(@(t,y) pendulum_ode(t,y,[ui;ui],[0;Ts],m,l,g,b,noise),[0,Ts],xi);
    sol = sol';

    % start from the next state
    xi = sol(:,end);

    % save trajectory
    state_trajectory = [state_trajectory,sol];
    time_trajectory = [time_trajectory,(i-1)*Ts+t'];
    control_trajory = [control_trajory,ui];
end

%% plot trajectory
figure;
tiledlayout(2,1)
nexttile
t_grid = 0:Ts:((num_steps-1)*Ts);
u_grid = control_trajory;
control_trajory_tt = interp1(t_grid,u_grid,time_trajectory,'previous');
plot(time_trajectory,control_trajory_tt,'LineWidth',2);
xlabel('$t$','FontSize',24,'Interpreter','latex');
ylabel('$u(t)$','FontSize',24,'Interpreter','latex');
ax = gca; ax.FontSize = 20;
grid on;

nexttile
plot(time_trajectory,state_trajectory,'LineWidth',2);
xlabel('$t$','FontSize',24,'Interpreter','latex');
ylabel('$x(t)$','FontSize',24,'Interpreter','latex');
legend('$\theta(t)$','$\dot{\theta}(t)$','FontSize',20,'Interpreter','latex')
ax = gca; ax.FontSize = 20;
grid on;
