% Code for Validating a Control Lyaponov Function 
% for a Pendulum with an Observer

clear; clc; close all; format compact;
syms b m l g real % damping, mass, length, gravity
syms theta theta_dot real % theta = output
syms theta_hat theta_dot_hat real
syms u real % u = control input
syms k1 k2 real % k1, k2 = observer gains
A = [0,1; 0, -b/(m*l^2)];
B = [0; (u - m*g*l*sin(theta))/(m*l^2)];
C = [1,0];
D = [C,0];
K = [k1; k2];

F = [A,zeros(2,2);zeros(2,2),(A - K*C)];
G = [B;0;0];
z0 = [pi;0;0;0];
x = [theta; theta_dot];
x_hat = [theta_hat; theta_dot_hat];
e = x_hat - x;
z = [x;e];
V = z'*z;
gradient(V);
Lf = dot(gradient(V),F*z);
Lg = dot(gradient(V),G);
