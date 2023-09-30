function dx = pendulum_ode(t,states,u_grid,t_grid,m,l,g,b,noise)
u_t = interp1(t_grid,u_grid,t); % piece-wise linear
x = states;
dx = [x(2);
      -1/(m*l^2) * (-u_t + b*x(2)+m*g*l*sin(x(1))) ];
% dx = dx + noise*randn(2,1);
dx = dx + [0; noise];
end