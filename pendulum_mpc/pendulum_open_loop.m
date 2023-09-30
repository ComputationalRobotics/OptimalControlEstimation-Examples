clc; clear; close all;

% pendulum parameters
m = 1; l = 1; g = 9.8; b = 0.1; umax = g/2;

initial_state = [0;0];
N = 101;
h = 0.1;

% trajectory optimization
[uopt,xopt] = PendulumTrajOpt(N,h,initial_state,m,l,g,b,umax,[],[]);

%% plot optimized plans
T = (N-1)*h;
t_grid = linspace(0,1,N)*T;

figure;
tiledlayout(2,1)
nexttile
scatter(t_grid,uopt,100,'filled'); hold on;
plot(t_grid,uopt,'LineWidth',2);
xlabel('$t$','FontSize',24,'Interpreter','latex');
ylabel('$u(t)$','FontSize',24,'Interpreter','latex');
ax = gca; ax.FontSize = 20;
grid on;

nexttile
scatter(t_grid,xopt,100,[0,0,1;1,0,0],'filled'); hold on;
plot(t_grid,xopt,'LineWidth',2);
xlabel('$t$','FontSize',24,'Interpreter','latex');
ylabel('$x(t)$','FontSize',24,'Interpreter','latex');
grid on;
legend('$\theta(t)$','$\dot{\theta}(t)$','FontSize',20,'Interpreter','latex');
ax = gca; ax.FontSize = 20;

%% deploy the controls
x = initial_state;
state_trajectory = [];
time_trajectory = [];
noise = 1.0;
for k = 1:N-1
    uk = uopt(k);
    [t,sol] = ode89(@(t,y) pendulum_ode(t,y,[uk;uk],[0;h],m,l,g,b,noise),[0,h],x);
    sol = sol';
    x = sol(:,end);
    % save trajectory
    state_trajectory = [state_trajectory,sol];
    time_trajectory = [time_trajectory,(k-1)*h+t'];
end

figure;
tiledlayout(1,1)
nexttile
plot(time_trajectory,state_trajectory,'LineWidth',2);
xlabel('$t$','FontSize',24,'Interpreter','latex');
ylabel('$x(t)$','FontSize',24,'Interpreter','latex');
grid on;
legend('$\theta(t)$','$\dot{\theta}(t)$','FontSize',20,'Interpreter','latex');
ax = gca; ax.FontSize = 20;


