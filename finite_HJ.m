%% parameters
dim = 1;
P = eye(dim);
Q = eye(dim);
R = 1;
T = 0.5;
x_max = 8;
N_mesh_x = 1000;
N_mesh_t = 1000;
A = 1;
B = 1;
%% solve it
% generate a mesh
mesh_x = linspace(-1,1,N_mesh_x)*x_max;
mesh_T = linspace(1,0,N_mesh_t)*T;
V = {};
% figure;
% terminal cost
V{end+1} = mesh_x.^2*P;
for i = 1:N_mesh_t-1
    V_cur = V{end};
    % gradient on x
    dV = zeros(1,N_mesh_x);
    dV(1) = (V_cur(2)-V_cur(1))/(mesh_x(2)-mesh_x(1));
    dV(end) = (V_cur(end)-V_cur(end-1))/(mesh_x(end)-mesh_x(end-1));
    for j = 2:N_mesh_x-1
        dV(j) = (V_cur(j+1)-V_cur(j-1))/(mesh_x(j+1)-mesh_x(j-1));
    end
    mesh_u = -B/2/R*dV;
    f = A*mesh_x + B*mesh_u;
    %gradient on t
    dV_t = -(mesh_x.^2*Q+mesh_u.^2*R+dV.*f);
    %simply Euler forward
    V_next = V_cur + dV_t*(mesh_T(i+1)-mesh_T(i));
    V{end+1} = V_next;
    % plot(mesh_x,V_next);
    % hold on
end
[X,Y] = meshgrid(mesh_x,mesh_T);
Z = zeros(N_mesh_t,N_mesh_x);

for i = 1:N_mesh_t
    Z(i,:) = V{i};
end
figure;
ax = gca; ax.FontSize = 20;
surf(X,Y,Z);
shading interp
xlabel('x','FontSize',24,'Interpreter','latex');
ylabel('t','FontSize',24,'Interpreter','latex');
ax = gca; ax.FontSize = 20;
%V_t = -(x^TQx + u^TRu + V_x^Tf(x,u))
%f(x,u) = Ax+Bu;
%% useful function
