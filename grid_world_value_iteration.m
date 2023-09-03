clc; clear; close all;

%% create problem data
% states
[X,Y] = meshgrid(1:10,1:10);
X = X(:);
Y = Y(:);
XY = [X,Y];
nx = length(X); % number of states

% controls
U = [1, 0;
     -1, 0;
     0, 1;
     0, -1;
     0, 0];
nu = size(U,1); % number of controls

% create vector g
target = [1,10]; 
obstacles = [2,2;
             2,3;
             3,2;
             3,3;
             3,4;
             4,5;
             4,6;
             5,5;
             5,6;
             6,5;
             6,6;
             7,7;
             7,8;
             8,2;
             8,3;
             8,4;
             8,7;
             8,8;
             9,7;
             9,8;
             10,7;
             10,8];

g = zeros(nx,nu);
for i = 1:nx
    for j = 1:nu
        x = XY(i,:); % state
        if ismember(x,target,'rows')
            state_cost = 0;
        elseif ismember(x,obstacles,'rows')
            state_cost = 20;
        else
            state_cost = 1;
        end

        u = U(j,:); % control
        control_cost = sum(u.^2);

        g(i,j) = state_cost;
    end
end
g_mat = g;
g = g(:);

% construct transition matrix P
P = zeros(nx*nu,nx);
count = 1;
for j = 1:nu
    u = U(j,:); % control
    for i = 1:nx
        x = XY(i,:); % state
        xp = x + u; % next state

        [~,loc] = ismember(xp,XY,'rows');

        if loc == 0 % if the next state is not in the grid, assume the state stays where it was
            p = zeros(1,nx);
            p(i) = 1;
            P(count,:) = p;
        else
            p = zeros(1,nx);
            p(loc) = 1;
            P(count,:) = p;
        end
        count = count + 1;
    end
end

%% start value iteration
Q = zeros(nx*nu,1);
iter = 1;
MAX_ITERS = 1e3;
while iter < MAX_ITERS
    Q_mat = reshape(Q,nx,nu);
    J_Q = min(Q_mat,[],2);

    Q_new = g + P*J_Q;
    Q_diff = norm(Q_new - Q);
    Q = Q_new;
    iter = iter + 1;

    % check convergence
    if Q_diff < 1e-6
        fprintf("Value iteration converged in %d iterations, Q_diff = %3.2e.\n",iter,Q_diff)
        break
    end
end
% compute optimal cost to go
Q_mat = reshape(Q,nx,nu);
J = min(Q_mat,[],2);

% plot optimal cost to go
J_mat = zeros(10,10);
for cnt = 1:nx
    x = XY(cnt,:);
    J_mat(x(1),x(2)) = J(cnt);
end
figure; image(J_mat,'CDataMapping','scaled'); colorbar;
set(gca,'XAxisLocation','top')
hold on

%% plot the optimal trajectory starting from a point
x = [8,5];
x_traj = x;
while true
    [~,loc] = ismember(x,XY,'rows');
    Q_x = Q_mat(loc,:);
    [~,id] = min(Q_x);
    u_star = U(id,:);
    x = x + u_star;
    x_traj = [x_traj;x];

    if norm(x - target) == 0
        break;
    end
end

plot(x_traj(:,2),x_traj(:,1),'r-x','LineWidth',2,'MarkerSize',10)
axis square
        

