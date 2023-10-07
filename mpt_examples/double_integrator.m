clc; clear; close all;
% requires the Multi-Parametric Toolbox (MPT) version 3
% https://www.mpt3.org/

% create double integrator system
A = [1, 1; 0, 1]; B = [0; 1];
sys = LTISystem('A',A,'B',B);
% add simple control and state bounds
sys.x.min = [-5;-5];
sys.x.max = [5;5];
sys.u.min = -0.5;
sys.u.max = 0.5;

% Compute the reachable set from the origin
% polyhedron of a single point
X0 = Polyhedron([0;0]);
S = sys.reachableSet('N',1);
S.plot()


