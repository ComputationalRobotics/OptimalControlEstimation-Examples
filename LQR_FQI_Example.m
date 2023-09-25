clc; clear; close all;

%% system parameters
h = 0.01;
A = [1, h; 0, 1];
B = [0; h];
Q = 0.1*eye(2);
R = 1;

%% groundtruth S matrix from LQR
[~,S_lqr,~] = dlqr(A,B,Q,R,zeros(2,1));
M_lqr = [Q+A'*S_lqr*A, A'*S_lqr*B;
         B'*S_lqr*A, B'*S_lqr*B + R];

%% Fitted Q-value iteration
num_samples     = 6;
num_iterations  = 1e4;
th              = 1e-8;
M_vec           = randn(6,1);
M_traj          = M_vec;
M_mat           = vec2mat(M_vec);
iter = 1;

while iter < num_iterations
    X = [];  % Feature matrix
    Y = [];  % Target cost vector
    
    for i = 1:num_samples
        x = randn(2,1);
        u = randn(1);

        ai = [A*x + B*u;0];
        L  = [0;0;1];
        
        % solve u'
        up = - ( (L'*M_mat*L) \ (L'*M_mat*ai)  );

        % compute beta
        beta = running_cost(x,u,Q,R) + (L*up+ai)'*M_mat*(L*up+ai);

        % Store data
        X = [X; [x(1)^2, 2*x(1)*x(2), 2*x(1)*u, x(2)^2, 2*x(2)*u, u^2]];
        Y = [Y; beta];
    end
    
    % Update weights
    M_vec_new = X \ Y;
    
    % Convergence check
    if norm(M_vec_new - M_vec) < th
        fprintf("FQI converged in %d iterations.\n",iter);
        break;
    end
    
    M_vec = M_vec_new;
    M_mat = vec2mat(M_vec);
    M_traj = [M_traj,M_vec];
    disp(M_vec);

    iter = iter + 1;
end

M_fqi = M_mat;

%% plot
labelsize = 20;
figure;
tiledlayout(3,1)
nexttile
plot(M_traj(1,:),'LineWidth',2)
ylabel('$M_1$','Interpreter','latex','FontSize',labelsize)
xlabel('Iteration','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;

nexttile
plot(M_traj(2,:),'LineWidth',2)
ylabel('$M_2$','Interpreter','latex','FontSize',labelsize)
xlabel('Iteration','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;

nexttile
plot(M_traj(3,:),'LineWidth',2)
ylabel('$M_3$','Interpreter','latex','FontSize',labelsize)
xlabel('Iteration','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;


%% helper functions
function g = running_cost(x,u,Q,R)
g = x'*Q*x + u'*R*u;
end

function M_mat = vec2mat(M_vec)
M_mat = [M_vec(1), M_vec(2), M_vec(3);
         M_vec(2), M_vec(4), M_vec(5);
         M_vec(3), M_vec(5), M_vec(6)];
end
