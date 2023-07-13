clc; clear; close all;

% Load solution from SOS program
sol = load('SOS-sols/cartpole_sol.mat').sol;

sostoolspath = "../SOSTOOLS";
addpath(genpath(sostoolspath))

m_c = 1; m_p = 1; l = 1; g = 9.8;
ns = 6;  % number of states
nu = 1;  % number of inputs
ny = 4;  % number of outputs

eps = 10^-2;

num_steps = 59990;
dt = 0.001;

%% simulate real dynamics

% State defined as x = [x, x_dot, sin(theta), cos(theta), theta_dot, a]
% where a = 1/(m_c/m_p+sin(theta)^2)

theta0 = pi/3;
s0 = sin(theta0);
c0 = cos(theta0);
a0 = 1/(m_c/m_p + sin(theta0^2));
x0 = [1;0;s0;c0;0;a0];
u_traj = zeros(nu,num_steps);
x_traj = zeros(ns,num_steps+1);
x_traj(:,1) = x0;
x = x0;
for i = 1:num_steps
    u = randi([-5, 5]);
    % u = 0;
    u_traj(i) = u;
    y = x([1,3,4,6]);
    xdot = cartpole_f(x,l) + cartpole_psi(u, y, m_p, l, m_c, g);

    x = x + xdot * dt;
    x(3:4) = normc(x(3:4));
    % TODO: check if a variation over time makes sense or not
    x_traj(:,i+1) = x;
end
y_traj = x_traj([1,3,4,6],:);
theta_traj = atan2(x_traj(3,:),x_traj(4,:));


%% simulate observer
xhat0 = [0;0;0;1;0;m_p/m_c];
xhat_traj = zeros(ns,num_steps+1);
xhat_traj(:,1) = xhat0;
xhat = xhat0;
for i = 1:num_steps
    yi = y_traj(:,i);
    u = u_traj(i);
    yhati = xhat([1,3,4,6]);
    Ce = yhati - yi;
    Qi = Q_func(Ce, eps, sol, ns);
    Mi = M_func(Ce, yi, sol, ns);
    
    xhatdot = cartpole_f(xhat,l) + cartpole_psi(u, y, m_p, l, m_c, g) + (Qi \ Mi)*Ce;

    xhat = xhat + xhatdot * dt;
    xhat(3:4) = normc(xhat(3:4));
    xhat_traj(:,i+1) = xhat;
end
theta_hat_traj = atan2(xhat_traj(3,:),xhat_traj(4,:));


%% plot comparison
labelsize = 20;
e_norm_traj = sqrt(sum((xhat_traj - x_traj).^2,1));
t_traj = (1:(num_steps+1)) * dt;

figure;

tiledlayout(5,1)
nexttile
x_comp = [x_traj(1,:);xhat_traj(1,:)];
plot(t_traj,x_comp','LineWidth',2)
ylabel('$x$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

nexttile
s_comp = [x_traj(3,:);xhat_traj(3,:)];
plot(t_traj,s_comp','LineWidth',2)
ylabel('$\sin \theta$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

nexttile
c_comp = [x_traj(4,:);xhat_traj(4,:)];
plot(t_traj,c_comp','LineWidth',2)
ylabel('$\cos \theta$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

nexttile
thdot_comp = [x_traj(5,:);xhat_traj(5,:)];
plot(t_traj,thdot_comp','LineWidth',2)
ylabel('$\dot{\theta}$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

nexttile
plot(t_traj,e_norm_traj','LineWidth',2)
ylabel('$\Vert \hat{x} - x \Vert$','Interpreter','latex','FontSize',labelsize)
xlabel('time','FontSize',labelsize)
ax = gca;
ax.FontSize = 16;

%% Create animation

filename = 'Animations/cartpole_animation.gif';
figure;
for i = 1:30:num_steps
    clf;
    
    % Calculate cart and pole positions
    cart_pos = x_traj(1, i);
    pole_pos = l * [sin(theta_traj(i)); -cos(theta_traj(i))];  % Pole is attached to the cart at one end and hangs down
    
    % Plot the cart
    rectangle('Position', [cart_pos-0.5, -0.25, 1, 0.5], 'FaceColor', 'b');
    hold on;
    
    % Plot the pole
    line([cart_pos, cart_pos+pole_pos(1)], [0, pole_pos(2)], 'Color', 'r', 'LineWidth', 2);
    
    xlims = [-2, 5];
    xlim(xlims);
    ylims = [-2, 2];
    ylim(ylims);
    axis equal;
    set(gca, 'YLimMode', 'manual', 'XLimMode', 'manual');
    xlabel('x');
    ylabel('y');

    % Add time annotation
    time = i * dt;  % calculate current time
    annotation_text = sprintf('Time: %.2f s', time);
    text(xlims(2)-1, ylims(1)+0.2, annotation_text, 'FontSize', 10);
    
    % Capture the plot as an image
    frame = getframe;
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF file
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', dt);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', dt);
    end
end


%% helper functions
function Q = Q_func(Ce, eps, sol, ns)
e = mpvar('e',[ns,1]);
e_1 = Ce(1);
e_3 = Ce(2);
e_4 = Ce(3);
e_6 = Ce(4);
Q_tilde = double(subs(sol.Q, e([1,3,4,6]), [e_1;e_3;e_4;e_6]));
Q = Q_tilde + eps*eye(ns);
end

function M = M_func(Ce, y, sol, ns)
e = mpvar('e',[ns,1]);
x = mpvar('x',[ns,1]);
e_1 = Ce(1);
e_3 = Ce(2);
e_4 = Ce(3);
e_6 = Ce(4);
M = double(subs(sol.M, [e([1,3,4,6]); x([1,3,4,6])], [e_1;e_3;e_4;e_6;y]));
end

% State defined as x = [x, x_dot, sin(theta), cos(theta), theta_dot, a]
% where a = 1/(m_c/m_p+sin(theta)^2)
function f = cartpole_f(x,l)
f = [x(2);
    x(6) * x(3)*l*x(5)^2;
    x(4)*x(5);
    -x(3)*x(5);
    -x(6)*x(5)^2*x(3)*x(4);
    -2*x(3)*x(4)*x(5)*x(6)^2;
    ];
end

function psi = cartpole_psi(u, y, m_p, l, m_c, g)
psi = [0;
    1/m_p*y(4)*u + y(4)*y(2)*y(3)*g;
    0;
    0;
    -1/(l*m_p)*y(4)*(u*y(3) + (m_c+m_p)*g*y(2));
    0;
    ];
end
