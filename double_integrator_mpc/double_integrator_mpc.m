clc; clear; close all;

cvxpath = "../../cvx";
addpath(genpath(cvxpath));


P = eye(2); Q = eye(2); R = 10; N = 3;
x1 = [-4.5;2];
x2 = [-4.5;3];

[u1_traj, x1_traj] = double_integrator_rhc(x1,P,Q,R,N);
[u2_traj, x2_traj] = double_integrator_rhc(x2,P,Q,R,N);

%% plot the results
figure;
scatter(x1_traj(1,:),x1_traj(2,:),100,'blue','filled');
hold on
plot(x1_traj(1,:),x1_traj(2,:),'blue','LineWidth',2);

scatter(x2_traj(1,:),x2_traj(2,:),100,'red','filled');
hold on
plot(x2_traj(1,:),x2_traj(2,:),'red','LineWidth',2);
grid on
xlabel('$x_1$','FontSize',24,'Interpreter','latex');
ylabel('$x_2$','FontSize',24,'Interpreter','latex');
ax = gca; ax.FontSize = 20;


%% Receding horizon control
function [u_traj,x_traj] = double_integrator_rhc(x,P,Q,R,N)
x_traj = x;
u_traj = [];
while true
    % solve open-loop convex optimization
    uopt = double_integrator_ocp(x,N,P,Q,R);
    % check feasibility
    if isempty(uopt)
        fprintf("Subproblem infeasible.\n")
        break;
    end
    u = uopt(1);
    u_traj = [u_traj;u];
    
    % go to next step
    x = double_integrator_dynamics(x,u);
    x_traj = [x_traj,x];

    % check convergence
    dist_to_origin = norm(x);
    if dist_to_origin < 1e-6
        fprintf("Current state distance to origin is %3.2e < %3.2e.\n",dist_to_origin);
        break
    end
end

end

