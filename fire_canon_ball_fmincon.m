clc; clear; close all;

g = 9.8;
x0 = [0;0;1]; % initial guess of the solution

% define objective function
obj = @(x) objective(x);

% define nonlinear constraints
nonlincon = @(x) nonlinear_con(x,g);

% define options for fmincon
options = optimoptions('fmincon','Algorithm','interior-point',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
    'checkGradients',false);

% call fmincon
[xopt,fopt,~,out] = fmincon(obj,x0,... % objective and initial guess
    -eye(3),zeros(3,1),... % linear inequality constraints
    [],[],... % no linear equality constraints
    [],[],... % no upper and lower bounds
    nonlincon,... % nonlinear constraints
    options);

fprintf("Maximum constraint violation: %3.2f.\n",out.constrviolation);
fprintf("Objective: %3.2f.\n",fopt);

% plot the solution
T = xopt(3);
t = 0:0.01:T;
xt = xopt(1)*t;
yt = xopt(2)*t - 0.5*g*(t.^2);
figure;
plot(xt,yt,'LineWidth',2);
hold on
scatter(10,10,100,"red",'filled','diamond');
axis equal; grid on;
xlabel('$x_1(t)$','FontSize',24,'Interpreter','latex');
ylabel('$x_2(t)$','FontSize',24,'Interpreter','latex');
ax = gca; ax.FontSize = 20;


%% helper functions
function [f,df] = objective(x)
% x is our decision variable: [v1, v2, T]

% objective function
f = 0.5 * (x(1)^2 + x(2)^2);

% gradient of the objective function (optional)
df = [x(1); x(2); 0]; % column vector
end

function [c,ceq,dc,dceq] = nonlinear_con(x,g)
% no inequality constraints
c = []; dc = [];

% two equality constraints
ceq = [x(1)*x(3)-10;
       x(2)*x(3)-0.5*g*x(3)^2-10];

% explicit gradient of ceq
dceq = [x(3), 0;
        0, x(3);
        x(1), x(2)-g*x(3)];

end
