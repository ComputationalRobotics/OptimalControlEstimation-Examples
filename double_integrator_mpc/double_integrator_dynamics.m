function xnew = double_integrator_dynamics(x,u)
xnew = [1, 1; 0, 1]*x + [0; 1]*u;
end