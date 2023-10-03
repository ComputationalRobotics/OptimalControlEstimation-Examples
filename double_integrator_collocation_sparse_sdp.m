clc; clear; close all;

spotpath = '../spotless';
sospath  = '../SOSprograms';
mosekpath= '/Users/hengyang/mosek';

addpath(genpath(spotpath));
addpath(genpath(sospath));
addpath(genpath(mosekpath));

initial_state = [-10;0];
final_state = [0;0];

N = 6;
d = 2*(N-2) + N + 1;
h = msspoly('h',1);
u = msspoly('u',N);
x = msspoly('x',2*(N-2));
v = [h;u;x];
x = reshape(x,2,N-2); % unknown states

x_full = [initial_state,x,final_state];

objective = (N-1)*h;
eqs = [];
for k=1:N-1
    uk = u(k);
    ukp1 = u(k+1);
    xk = x_full(:,k);
    xkp1 = x_full(:,k+1);
    fk = double_integrator(xk,uk);
    fkp1 = double_integrator(xkp1,ukp1);
    
    % collocation points
    xkc = 0.5*(xk+xkp1) + h/8 * (fk - fkp1);
    ukc = 0.5*(uk + ukp1);
%     dxkc = -3/(2*h) * (xk-xkp1) - 0.25*(fk + fkp1);
    dxkc_h = -3/2 * (xk-xkp1) - 0.25*h*(fk + fkp1);
    
    % collocation constraint
    eqs = [eqs;
           dxkc_h - h*double_integrator(xkc,ukc)];
end

% array of cliques
cliques = {[1,2,3,N+1+blkIndices(1,2)]};
for k=2:N-2
    cliques = [cliques;
               {[1,k+1,k+2,N+1+blkIndices(k-1,2),N+1+blkIndices(k,2)]}];

end
cliques = [cliques;
           {[1,N,N+1,N+1+blkIndices(N-2,2)]}];

% array of equality constraint indices assigned to each clique
I = {};
for k = 1:N-1
    I = [I; 
         {blkIndices(k,2)}];
end

% array of inequality constraint indices assigned to each clique
ineqs = [h];
for k=1:N
    ineqs = [ineqs;
             1 - u(k)^2];
end
for i = 1:length(cliques)
    ids = cliques{i};
    vars = v(ids);
    ineqs = [ineqs;
             200 - vars'*vars];
end
J = {};
for i = 1:length(cliques)
    J{i} = [1,i+1,i+2,N+1+i];
end

orders = 2*ones(1,length(cliques));
[blk, At, C, b] = sparsesos(-1, [objective;-1], ...
    ineqs, eqs, v, ...
    cliques, J, I, orders);

% [At,b,c,k] = SDPT3data_SEDUMIdata(blk,At,C,b);
prob       = convert_sedumi2mosek(SDPT3data_SEDUMIdata(blk,At,C,b));
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);

figure; bar(eig(Xopt{1}));


%% helper functions
function dx = double_integrator(x,u)
A = [0 1; 0 0];
B = [0; 1];
dx = A * x + B * u;
end