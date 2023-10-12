clc; clear; close all;

A = [1,1;0,1]; B = [0;1];
Q = [1,0;0,0]; R = 1; 

[~,S,~] = dlqr(A,B,Q,R,zeros(2,1));

P = S;
N = 1;

x = [5;10];
x_traj = x;
u_traj = [];
for t = 1:100
    [uopt,info] = double_integrator_rhc(A,B,P,Q,R,N,x);

    xnew = A*x + B*uopt(1);

    if norm(xnew) < 1e-6
        fprintf("Close to origin.\n")
        break;
    end

    if norm(xnew - x) < 1e-9
        fprintf("Slow convergence.\n");
        break;
    end
    
    x = xnew;
    x_traj = [x_traj,x];
    u_traj = [u_traj,uopt(1)];
end

figure; plot(x_traj');

function [uopt,info] = double_integrator_rhc(A,B,P,Q,R,N,xt)
% decision variable
% [u(0),...,u(N-1),x(1),...,x(N)]
lb = [-1*ones(N,1);-inf*ones(2*N,1)];
ub = [1*ones(N,1);inf*ones(2*N,1)];

Aeq = [];
beq = [];
for i = 1:N-1
    Ai = zeros(2,3*N);
    Ai(:,i+1) = B;
    Ai(:,N+blkIndices(i,2)) = A;
    Ai(:,N+blkIndices(i+1,2)) = -eye(2);

    bi = zeros(2,1);

    Aeq = [Aeq;Ai];
    beq = [beq;bi];
end

A0 = zeros(2,3*N);
A0(:,1) = -B;
A0(:,N+blkIndices(1,2)) = eye(2);
b0 = A*xt;

Aeq = [A0;Aeq]; beq = [b0;beq];
blks = {};
for i = 1:N
    blks = [blks;{R}];
end
for i = 1:N-1
    blks = [blks;{Q}];
end
blks = [blks;{P}];

H = blkdiag(blks{:});

xopt = quadprog(H,[],[],[],Aeq,beq,lb,ub);

uopt = xopt(1:N);
info.H = H;
info.Aeq = Aeq;
info.beq = beq;
info.xopt = xopt;
end



