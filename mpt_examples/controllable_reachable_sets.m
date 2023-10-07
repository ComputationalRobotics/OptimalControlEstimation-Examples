clc; clear; close all;
% requires the Multi-Parametric Toolbox (MPT) version 3
% https://www.mpt3.org/
% Installation: https://www.mpt3.org/Main/Installation
% Example 10.3 and 10.4 in BBM book

% define LTI system
A = [1.5, 0; 1, -1.5]; B = [1; 0];
sys = LTISystem('A',A,'B',B);

% define state constraint
calX = Polyhedron('A',...
    [1,0;0,1;-1,0;0,-1], ...
    'b',[10;10;10;10]);
% define control constraint
calU = Polyhedron('A',[1;-1],'b',[5;5]);

%% compute controllable sets
% target set
S = Polyhedron('A',[1,0;0,1;-1,0;0,-1],'b',[1;1;1;1]);
K = [S];
N = 4;
for i = 1:N
    Si = sys.reachableSet('X',K(i),'U',calU,'N',1,'direction','backwards');
    K = [K; intersect(Si,calX)];
end

% plot (unshifted) controllable sets
figure; K.plot('Alpha',0.4);
legend('$\mathcal{K}_0(\mathcal{S})$',...
    '$\mathcal{K}_1(\mathcal{S})$',...
    '$\mathcal{K}_2(\mathcal{S})$',...
    '$\mathcal{K}_3(\mathcal{S})$',...
    '$\mathcal{K}_4(\mathcal{S})$','FontSize',20,'Interpreter','latex');
ax = gca; ax.FontSize = 16;
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_2$','FontSize',20,'Interpreter','latex');
% plot shifted controllable sets
Kshift = [];
for i = 1:length(K)
    Kshift = [Kshift; K(i) + (i-1)*[-8;0]];
end
figure; Kshift.plot();
legend('$\mathcal{K}_0(\mathcal{S})$',...
    '$\mathcal{K}_1(\mathcal{S})$',...
    '$\mathcal{K}_2(\mathcal{S})$',...
    '$\mathcal{K}_3(\mathcal{S})$',...
    '$\mathcal{K}_4(\mathcal{S})$','FontSize',20,'Interpreter','latex');
ax = gca; ax.FontSize = 16;
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_2$','FontSize',20,'Interpreter','latex');

%% Compute reachable sets
% initial set
X0 = Polyhedron('A',[1,0;0,1;-1,0;0,-1],'b',[1;1;1;1]);
R = [X0];
N = 4;
for i = 1:N
    Ri = sys.reachableSet('X',R(i),'U',calU,'N',1,'direction','forward');
    R = [R; intersect(Ri,calX)];
end
% plot (unshifted) reachable sets
figure; R.plot('Alpha',0.4);
legend('$\mathcal{R}_0(\mathcal{X}_0)$',...
    '$\mathcal{R}_1(\mathcal{X}_0)$',...
    '$\mathcal{R}_2(\mathcal{X}_0)$',...
    '$\mathcal{R}_3(\mathcal{X}_0)$',...
    '$\mathcal{R}_4(\mathcal{X}_0)$','FontSize',20,'Interpreter','latex');
ax = gca; ax.FontSize = 16;
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_2$','FontSize',20,'Interpreter','latex');
% plot shifted reachable sets
Rshift = [];
for i = 1:length(R)
    Rshift = [Rshift; R(i) + (i-1)*[40;0]];
end
figure; Rshift.plot();
legend('$\mathcal{R}_0(\mathcal{X}_0)$',...
    '$\mathcal{R}_1(\mathcal{X}_0)$',...
    '$\mathcal{R}_2(\mathcal{X}_0)$',...
    '$\mathcal{R}_3(\mathcal{X}_0)$',...
    '$\mathcal{R}_4(\mathcal{X}_0)$','FontSize',20,'Interpreter','latex');
ax = gca; ax.FontSize = 16;
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_2$','FontSize',20,'Interpreter','latex');