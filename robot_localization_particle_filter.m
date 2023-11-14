clc; clear; close all;

N = 20; % number of landmarks
landmarks = randn(2,N); % random landmarks
sig_w = 0.01;
sig_v = 0.01;
M = 1000; % number of samples
dt = 0.1; % sampling time


s0 = [0,0,0]; % initial state of the vehicle
calX0 = mvnrnd(ones(3,1),0.1*eye(3),M)'; % initial set of particles

% number of steps
num_steps = 200;

s = s0;
calX = calX0;
s_traj = zeros(3,num_steps);
s_traj(:,1) = s0;
calX_traj = zeros(3,M,num_steps);
calX_traj(:,:,1) = calX0;
alpha = zeros(M,1);

figure;

for i = 1:num_steps
    % generate true next state
    a = randn(2,1);
    s = vehicle(s,a,dt) + sig_w * randn(3,1);
    % generate measurement
    o = measure(s,landmarks);

    %% particle filter
    % propagate particles
    calX_p = zeros(3,M);
    for j = 1:M
        calX_p(:,j) = vehicle(calX(:,j),a,dt) + sig_w * randn(3,1);
    end
    % compute measurement errors 
    errs = zeros(N,M);
    for j = 1:M
        errs(:,j) = o - measure(calX_p(:,j),landmarks);
    end
    % compute weights
    alpha = exp(-0.5*sum(errs.^2,1)); 
    alpha = alpha / (sum(alpha)); alpha = alpha(:);
    % importance sampling
    ids = randsample(M,M,true,alpha);
    calX = calX_p(:,ids);

    clf;
    scatter(landmarks(1,:),landmarks(2,:),150,"black",'filled','square');
    hold on
    scatter(calX(1,:),calX(2,:),30,'blue','filled','o');
    scatter(s(1),s(2),150,"red",'filled','diamond');
    xlim([-2,2])
    ylim([-2,2])
    axis equal
    mytitle = sprintf("Particle Filter: Step %d",i);
    title(mytitle,'FontSize',22)
    pause(0.1)

end
    


%% helper functions
function snew = vehicle(s,a,dt)
x = s(1); y = s(2); theta = s(3);
l = a(1); u = a(2);
xnew = x + l*cos(theta)*dt;
ynew = y + l*sin(theta)*dt;
thetanew = theta + u*dt;
snew = [xnew; ynew; thetanew];
end

function o = measure(s,landmarks)
N = size(landmarks,2);
o = zeros(N,1);
for i = 1:N
    o(i) = norm(landmarks(:,i) - s(1:2));
end
end