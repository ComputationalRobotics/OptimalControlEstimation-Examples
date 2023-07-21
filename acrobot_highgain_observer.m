clc; clear; close all;
% 
cvxpath = "./cvx";
addpath(genpath(cvxpath))

ns = 4; ny = 2; nu = 1;
m_1 = 1; m_2 = 1; l_1 = 1; l_2 = 1; l_c1 = l_1/2; l_c2 = l_2/2; g = 9.8;
I_1 = 1/3*m_1*l_1^2/3; I_2 = 1/3*m_2*l_2^2;

obs_idx = [1, 2];

%% Find the optimal K that maximizes convergence rate
% TODO
k_1 = 1; k_2 = 0.1; L = 1000;

%% simulate real dynamics
dt = 0.0001;
T = 10;
num_steps = floor(T/dt);

x0 = [pi/6; pi/6; 1; 1];
u_traj  = zeros(nu,num_steps);
x_traj = zeros(ns,num_steps+1);
x_traj(:,1) = x0;
x = x0;
for i = 1:num_steps
    u     = rand/3;
    u_traj(:,i) = u;
    y     = x(obs_idx);
    xdot1 = x(3:4);
    xdot2 = acrobot_Phi(x,u,m_1,m_2,I_1,I_2,l_1,l_c1,l_c2,g);
    xdot  = [xdot1; xdot2];
    x     = x + xdot * dt;
    x_traj(:,i+1) = x;
end
y_traj = x_traj(obs_idx,:);


%% simulate observer
xhat0 = zeros(ns,1);
xhat_traj = zeros(ns,num_steps+1);
xhat_traj(:,1) = xhat0;
xhat = xhat0;
for i = 1:num_steps
    y = y_traj(:,i);
    u = u_traj(:,i);
    xhatdot1 = xhat(3:4) - L*k_1*(xhat(1:2) - y);
    xhatdot2 = acrobot_Phi(xhat,u,m_1,m_2,I_1,I_2,l_1,l_c1,l_c2,g) - L^2*k_2*(xhat(1:2) - y);
    xhatdot  = [xhatdot1; xhatdot2];
    xhat = xhat + xhatdot * dt;
    xhat_traj(:,i+1) = xhat;
end
e_norm_traj = sqrt(sum((xhat_traj - x_traj).^2,1));


%% plot comparison
labelsize = 20;

PLOT = true;
if PLOT
    t_traj = (1:(num_steps+1)) * dt;
    figure('Position', [300 200 900 900]);
    tiledlayout(5,1)
    nexttile
    plot(t_traj,x_traj(1,:),'LineWidth',2)
    hold on
    plot(t_traj,xhat_traj(1,:),'LineWidth',2)
    ylabel('$\theta_1$','Interpreter','latex','FontSize',labelsize)
    xlabel('time','FontSize',labelsize)
    legend('$\theta_1$','$\hat{\theta}_1$','Interpreter','latex','FontSize',labelsize)
    ax = gca;
    ax.FontSize = 16;

    nexttile
    plot(t_traj,x_traj(2,:),'LineWidth',2)
    hold on
    plot(t_traj,xhat_traj(2,:),'LineWidth',2)
    ylabel('$\theta_2$','Interpreter','latex','FontSize',labelsize)
    xlabel('time','FontSize',labelsize)
    legend('$\theta_2$','$\hat{\theta}_2$','Interpreter','latex','FontSize',labelsize)
    ax = gca;
    ax.FontSize = 16;

    nexttile
    plot(t_traj,x_traj(3,:),'LineWidth',2)
    hold on
    plot(t_traj,xhat_traj(3,:),'LineWidth',2)
    ylabel('$\dot{\theta}_1$','Interpreter','latex','FontSize',labelsize)
    xlabel('time','FontSize',labelsize)
    legend('$\dot{\theta}_1$','$\dot{\hat{\theta}}_1$','Interpreter','latex','FontSize',labelsize)
    ax = gca;
    ax.FontSize = 16;

    nexttile
    plot(t_traj,x_traj(4,:),'LineWidth',2)
    hold on
    plot(t_traj,xhat_traj(4,:),'LineWidth',2)
    ylabel('$\dot{\theta}_2$','Interpreter','latex','FontSize',labelsize)
    xlabel('time','FontSize',labelsize)
    legend('$\dot{\theta}_2$','$\dot{\hat{\theta}}_2$','Interpreter','latex','FontSize',labelsize)
    ax = gca;
    ax.FontSize = 16;

    nexttile
    e_norm_traj_log = log(e_norm_traj);
    plot(t_traj,e_norm_traj_log,'LineWidth',2)
    ylabel('$\log \Vert \hat{x} - x \Vert$','Interpreter','latex','FontSize',labelsize)
    xlabel('time','FontSize',labelsize)
    ax = gca;
    ax.FontSize = 16;
end

%% Create animation

ANIMATE = false;
if ANIMATE
    filename = 'Animations/acrobot_animation.gif';
    figure;

    delay = round(10^5 * dt);
    for i = 1:delay:num_steps
        clf;
        
        % Calculate link positions
        link1_pos = [l_1 * sin(x_traj(1, i)), -l_1 * cos(x_traj(1, i))];
        link2_pos = link1_pos + [l_2 * sin(x_traj(1, i) + x_traj(2, i)), -l_2 * cos(x_traj(1, i) + x_traj(2, i))];
        
        % Plot the links
        hold on;
        plot([0, link1_pos(1)], [0, link1_pos(2)], 'b', 'LineWidth', 2);
        plot([link1_pos(1), link2_pos(1)], [link1_pos(2), link2_pos(2)], 'r', 'LineWidth', 2);
        plot(0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');  % Origin point
        xlims = [-l_1-l_2-1, l_1+l_2+1];
        ylims = [-l_1-l_2-1, l_1+l_2+1];
        xlim(xlims);
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
end


%% Helper functions
% State defined as x = [theta_1; theta_2; thetadot_1; thetadot_2]
function Phi = acrobot_Phi(x,u,m_1,m_2,I_1,I_2,l_1,l_c1,l_c2,g)
    M = [I_1 + I_2 + m_2*l_1^2 + 2*m_2*l_1*l_c2*cos(x(2)), I_2 + m_2*l_1*l_c2*cos(x(2));
         I_2 + m_2*l_1*l_c2*cos(x(2)), I_2];
    C = [-m_2*l_1*l_c2*sin(x(2))*x(4), -m_2*l_1*l_c2*sin(x(2))*x(4);
         m_2*l_1*l_c2*sin(x(2))*x(3), 0];
    g = [-m_1*l_c1*g*sin(x(1)) - m_2*g*(l_1*sin(x(1)) + l_c2*sin(x(1)+x(2)));
         -m_2*l_c2*g*sin(x(1)+x(2))];
    tau = [0; u];
    Phi = M \ (tau - C*x(3:4) + g);
end
