clc; clear; close all;
% 
cvxpath = "../cvx";
addpath(genpath(cvxpath))

m = 1; g = 9.8; l = 1; b = 0.1;
A = [0, 1; 0, -b/(m*l^2)];
C = [1, 0];

%% Estimate convergence rate gamma for varying k
k_choices = [0.1,1,10,100,1000,10000];
gams = zeros(length(k_choices),1);
for i = 1:length(k_choices)
    k = k_choices(i);
    K = [k; 0];
    AKC = A - K*C;
    cvx_begin
        variable P(2,2) symmetric 
        minimize(0)
        subject to
            P == semidefinite(2) 
            AKC'*P + P*AKC == - eye(2)
    cvx_end
    max_eig_P = max(eig(P));
    gams(i) = 0.5 / max_eig_P;
end

%% Find the optimal K that maximizes convergence rate
cvx_begin
    variable H(2)
    variable P(2,2) symmetric
    minimize(lambda_max(P))
    subject to
        P == semidefinite(2)
        A'*P - C'*H' + P*A - H*C == - eye(2);
cvx_end
max_gam = 0.5 / max(eig(P));
K = P \ H;