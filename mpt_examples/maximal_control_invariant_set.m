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

Omega = [calX];
while true
    last_Omega = Omega(end);
    pre_omega = sys.reachableSet('X',last_Omega,'U',calU,'N',1,...
        'direction','backwards');
    new_Omega = intersect(pre_omega,last_Omega);
    Omega = [Omega;new_Omega];
    if new_Omega == last_Omega
        fprintf("Converged to maximal control invariant set.\n");
        break;
    end
end
figure; Omega.plot('Color','gray','Alpha',0.2);
ax = gca; ax.FontSize = 16;
xlabel('$x_1$','FontSize',20,'Interpreter','latex');
ylabel('$x_2$','FontSize',20,'Interpreter','latex');

%% use MPT function
C = sys.invariantSet('X', calX, 'U', calU);
figure; C.plot();

    