clc; clear; close all;

%% system parameters
h = 0.01;
A = [1, h; 0, 1];
B = [0; h];
Q = 0.1*eye(2);
R = 1;

%% groundtruth S matrix from LQR
[~,S_lqr,~] = dlqr(A,B,Q,R,zeros(2,1));

%% Fitted value iteration
num_samples     = 2;
num_iterations  = 1e4;
th              = 1e-8;
S_vec           = randn(3,1);
S_traj          = S_vec;
S_mat           = vec2mat(S_vec);
iter = 1;

while iter < num_iterations
    X = [];  % Feature matrix
    Y = [];  % Target cost vector
    
    for i = 1:num_samples
        x = randn(2,1);
        
        % solve u
        u = - ( (1 + B'*S_mat*B) \ B'*S_mat*A*x );

        % compute beta
        beta = running_cost(x,u,Q,R) + (A*x+B*u)'*S_mat*(A*x+B*u);

        % Store data
        X = [X; [x(1)^2, 2*x(1)*x(2), x(2)^2]];
        Y = [Y; beta];
    end
    
    % Update weights
    S_vec_new = X \ Y;
    
    % Convergence check
    if norm(S_vec_new - S_vec) < th
        fprintf("FVI converged in %d iterations.\n",iter);
        break;
    end
    
    S_vec = S_vec_new;
    S_mat = vec2mat(S_vec);
    S_traj = [S_traj,S_vec];
    disp(S_vec);

    iter = iter + 1;
end

S_fvi = S_mat;

%% plot
labelsize = 20;
figure;
tiledlayout(3,1)
nexttile
plot(S_traj(1,:),'LineWidth',2)
ylabel('$S_1$','Interpreter','latex','FontSize',labelsize)
xlabel('Iteration','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;

nexttile
plot(S_traj(2,:),'LineWidth',2)
ylabel('$S_2$','Interpreter','latex','FontSize',labelsize)
xlabel('Iteration','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;

nexttile
plot(S_traj(3,:),'LineWidth',2)
ylabel('$S_3$','Interpreter','latex','FontSize',labelsize)
xlabel('Iteration','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;


%% helper functions
function g = running_cost(x,u,Q,R)
g = x'*Q*x + u'*R*u;
end

function S_mat = vec2mat(S_vec)
S_mat = [S_vec(1), S_vec(2);
         S_vec(2), S_vec(3)];
end
