clc; clear; close all;

m = 1; gravity = 9.8; l = 1; b = 0.1;

N = 50; % number of bins in discretization

THETA       = linspace(-pi,pi,N); % theta between -pi and pi
THETADOT    = linspace(-pi,pi,N); % thetadot between -pi and pi
U           = linspace(-4.9,4.9,N)'; % control between -4.9 and 4.9

[X1,X2] = meshgrid(THETA,THETADOT);
X = [X1(:),X2(:)];

nx = size(X,1); % number of states
nu = length(U); % number of controls

% create running cost vector
g = pendulum_running_cost(X,U);
g_mat = reshape(g,nx,nu);

% construct transition matrix P
% create a triangulation of the mesh points
DT = delaunayTriangulation(X);
figure; triplot(DT);
dt = 0.01; % convert continuous-time dynamics to discrete-time
P_rows = []; % row indices of nonzero entries in P
P_cols = []; % col indices of nonzero entries in P
P_vals = []; % values of nonzero entries in P
count = 1;
for j = 1:nu
    u = U(j); % control
    for i = 1:nx
        x = X(i,:); % state
        xp = pendulum_f_z(x(:),u,m,gravity,l,b)*dt + x(:); % next state
        % wrap theta between -pi and pi
        xp(1) = clip(xp(1),-pi,pi); % use wrapToPi also works
        % clip thetadot between -pi and pi
        xp(2) = clip(xp(2),-pi,pi);
        % find xp's barycentric coordinate in the triangulation
        [ID,B] = pointLocation(DT,xp');
        % create probability matrix
        if length(ID) > 1
            error("point lies inside two triangles.")
        end
        P_rows = [P_rows,count,count,count];
        P_cols = [P_cols,DT.ConnectivityList(ID,:)];
        P_vals = [P_vals,B];
        count = count + 1;
    end
end
P = sparse(P_rows,P_cols,P_vals,nx*nu,nx);

%% start value iteration
fprintf("Start value iteration ...")
gamma = 0.999;
Q = zeros(nx*nu,1);
iter = 1;
MAX_ITERS = 1e5;
while iter < MAX_ITERS
    Q_mat = reshape(Q,nx,nu);
    J_Q = min(Q_mat,[],2);

    Q_new = g + gamma*P*J_Q;
    Q_diff = norm(Q_new - Q);
    Q = Q_new;
    iter = iter + 1;
    % check convergence
    if Q_diff < 1e-9
        fprintf("Value iteration converged in %d iterations, Q_diff = %3.2e.\n",iter,Q_diff)
        break
    end
end

%% compute optimal cost to go
Q_mat = reshape(Q,nx,nu);
J = min(Q_mat,[],2);

figure;
surf(reshape(X(:,1),N,N),...
     reshape(X(:,2),N,N),...
     reshape(J,N,N));
xlabel("$\theta$",'Interpreter','latex','FontSize',16);
ylabel("$\dot{\theta}$",'Interpreter','latex','FontSize',16);
zlabel("$J^\star$",'Interpreter','latex','FontSize',16)

%% execute the controller
x = [-0.4;0.3];
x_traj = x;
u_traj = [];
MAX_ITERS = 1e3;
iter = 1;
while iter < MAX_ITERS
    [ID,B] = pointLocation(DT,x');
    vertices = DT.ConnectivityList(ID,:);
    u_v = zeros(length(vertices),1);
    for j = 1:length(vertices)
        v = vertices(j);
        [~,id] = min(Q_mat(v,:));
        u_v(j) = U(id);
    end
    u = B(:)'*u_v;
    xp = pendulum_f_z(x(:),u,m,gravity,l,b)*dt + x(:); % next state
    u_traj = [u_traj;u];

    if norm(x - xp) < 1e-6
        break;
    end    
    x = xp;
    x_traj = [x_traj,x];
    iter = iter + 1;
end

%% plot comparison
labelsize = 20;

t_traj = (1:(size(x_traj,2)-1)) * dt;
figure;
tiledlayout(3,1)
nexttile
plot(t_traj,x_traj(1,1:end-1),'LineWidth',2)
ylabel('$z_1$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;

nexttile
plot(t_traj,x_traj(2,1:end-1),'LineWidth',2)
ylabel('$z_2$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;

nexttile
plot(t_traj,u_traj(1:length(t_traj)),'LineWidth',2)
ylabel('$u$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
grid on
ax = gca;
ax.FontSize = 16;


%% helper functions
% pendulum continuous-time dynamics
function zdot = pendulum_f_z(z,u,m,g,l,b)
z1 = z(1);
z2 = z(2);

zdot = [
    z2;
    1/(m*l^2)*(u - b*z2 + m*g*l*sin(z1))
];
end

% create pendulum g vector
function g = pendulum_running_cost(X,U)
nx = size(X,1); 
nu = length(U);

g_x = sum(X.^2,2);
g_u = U.^2;

g = repmat(g_x,1,nu) + repmat(g_u',nx,1);

g = g(:);
end

function u_clip = clip(u,umin,umax)
u_clip = min(umax,max(u,umin));
end



