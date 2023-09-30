function [uopt,xopt] = PendulumTrajOpt(N,h,initial_state,m,l,g,b,umax,u_guess,x_guess)
if isempty(x_guess)
    x0 = randn(2,N);
else
    x0 = x_guess;
end

if isempty(u_guess)
    u0 = randn(N,1);
else
    u0 = u_guess;
end

v0 = [x0(:); u0(:)];
lb = [-Inf*ones(2*N,1);
      -umax*ones(N,1)];
ub = [Inf*ones(2*N,1);
      umax*ones(N,1)];

obj = @(v) objective(v,N,h);
nonlincon = @(v) collocation(v,N,h,initial_state,m,l,g,b);

options = optimoptions('fmincon','Algorithm','interior-point',...
    'display','iter','MaxFunctionEvaluations',1e6,'MaxIterations',1e4);

[vopt,fopt,~,out] = fmincon(obj,v0,...
    [],[],[],[],... % no linear (in)equality constraints
    lb,ub,...
    nonlincon,options);

fprintf("Maximum constraint violation: %3.2f.\n",out.constrviolation);
fprintf("Objective: %3.2f.\n",fopt);

uopt = vopt(2*N+1:end);
xopt = vopt(1:2*N);
xopt = reshape(xopt,2,N);

end


%% helper functions
function f = objective(v,N,h)
x = v(1:2*N); x = reshape(x,2,N);
theta = x(1,:);
thetadot = x(2,:);
u = v(2*N+1:end);

tmp = [cos(theta);
       sin(theta);
       thetadot];

x_cost = sum((tmp - [-1;0;0]).^2,1);
u_cost = u.^2;
f = sum(x_cost(:) + u_cost(:)) * h;
end


function dx = pendulum(x,u,m,l,g,b)
    dx = [x(2);
          -1/(m*l^2) * (-u + b*x(2)+m*g*l*sin(x(1))) ];
end

function [c,ceq] = collocation(v,N,h,initial_state,m,l,g,b)
x = reshape(v(1:2*N),2,N);
u = v(2*N+1:end);

c = [];
ceq = [];

for k=1:N-1
    uk = u(k);
    ukp1 = u(k+1);
    xk = x(:,k);
    xkp1 = x(:,k+1);
    fk = pendulum(xk,uk,m,l,g,b);
    fkp1 = pendulum(xkp1,ukp1,m,l,g,b);
    
    % collocation points
    xkc = 0.5*(xk+xkp1) + h/8 * (fk - fkp1);
    ukc = 0.5*(uk + ukp1);
    dxkc = -3/(2*h) * (xk-xkp1) - 0.25*(fk + fkp1);
    
    % collocation constraint
    ceq = [ceq;
           dxkc - pendulum(xkc,ukc,m,l,g,b)];
end

ceq = [ceq;
       x(:,1) - initial_state;
       x(:,end) - [pi;0]];
end


