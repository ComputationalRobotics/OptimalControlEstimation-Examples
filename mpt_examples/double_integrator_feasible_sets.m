clc; clear; close all;
% requires the Multi-Parametric Toolbox (MPT) version 3
% https://www.mpt3.org/

% create double integrator system
A = [1, 1; 0, 1]; B = [0; 1];
sys = LTISystem('A',A,'B',B);

% define state constraint
calX = Polyhedron('A',...
    [1,0;0,1;-1,0;0,-1], ...
    'b',[5;5;5;5]);
% define control constraint
calU = Polyhedron('A',[1;-1],'b',[0.5;0.5]);
N = 3;

%% compute feasible sets
% Xf = calX;
Xf = Polyhedron('A',...
    [1,0;0,1;-1,0;0,-1], ...
    'b',[0;0;0;0]);
F = [Xf];
for i = 1:N
    last_X = F(end);
    new_X = intersect(calX,...
        sys.reachableSet('X',last_X,'U',calU,'N',1,'direction','backwards'));
    F = [F;new_X];
end
colors = [0,0,0;1,0,0;0,1,0;0,0,1];
figure; F.plot('Alpha',0.4,'ColorMap',colors);
legend('$\mathcal{X}_3$',...
    '$\mathcal{X}_2$',...
    '$\mathcal{X}_1$',...
    '$\mathcal{X}_0$','FontSize',20,'Interpreter','latex');
ax = gca; ax.FontSize = 16;
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_2$','FontSize',20,'Interpreter','latex');


