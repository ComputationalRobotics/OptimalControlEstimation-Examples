clc
clear
close all

g = [1;0;1;1;1;0;1;1];
P = [1,0,0,0;
     0,1,0,0;
     1,0,0,0;
     0,1,0,0;
     0,1,0,0;
     0,1,0,0;
     0,0,0,1;
     0,0,0,1];
P = sparse(P);

J = rand(4,1);

MAX_ITERS = 1e3;
iter = 1;
while iter < MAX_ITERS
    tmp = P * J + g;
    tmp = reshape(tmp,4,2);
    Jnew = min(tmp,[],2);
    if norm(Jnew - J) < 1e-6
        break;
    end
    J = Jnew;
    iter = iter + 1;
end
