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
v = msspoly('v',d);
h = v(1); % time interval
u = v(2:1+N); % unknown controls
x = reshape(v(N+2:end),2,N-2); % unknown states

x_full = [initial_state,x,final_state];

problem.vars = v;
problem.objective = (N-1)*h;
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
problem.equality = eqs; 

ineqs = [h;100*(N+1)+N - v'*v];
for k=1:N
    ineqs = [ineqs;
             1 - u(k)^2];
end

problem.inequality = ineqs;
kappa              = 2; % relaxation order
[SDP,info]         = dense_sdp_relax(problem,kappa);

prob       = convert_sedumi2mosek(SDP.sedumi.At,...
                                  SDP.sedumi.b,...
                                  SDP.sedumi.c,...
                                  SDP.sedumi.K);
[~,res]    = mosekopt('minimize info',prob);
[Xopt,yopt,Sopt,obj] = recover_mosek_sol_blk(res,SDP.blk);

figure; bar(eig(Xopt{1}));


%% helper functions
function dx = double_integrator(x,u)
A = [0 1; 0 0];
B = [0; 1];
dx = A * x + B * u;
end