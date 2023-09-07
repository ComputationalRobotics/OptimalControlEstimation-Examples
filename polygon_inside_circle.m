clc; clear; close all;

% number of points to be placed
N = 3;

% define the objective function
% fmincon assumes minimization
% We minimize the negative perimeter so as to maximize the perimeter
objective = @(u) -1*perimeter(u);

% choose which algorithm to use for solving
options = optimoptions('fmincon', 'Algorithm', 'interior-point');

% supply an initial guess
% since this is a convex problem, we can use any initial guess
u0 = rand(N,1);

% solve
uopt = fmincon(objective,u0,... % objective and initial guess
    -eye(N),zeros(N,1),... % linear inequality constraints
    ones(1,N),2*pi,... % linear equality constraints
    [],[],[],... % we do not have lower/upper bounds and nonlinear constraints
    options);

% plot the solution
x = zeros(N,1);
for k = 2:N
    x(k) = x(k-1) + uopt(k-1);
end
figure;
% plot a circle
viscircles([0,0],1);
hold on
% scatter the placed points
scatter(cos(x),sin(x),'blue','filled');
axis equal;

%% helper functions
% The objective function
function f = perimeter(u)
f = sum(2*sin(u/2));
end