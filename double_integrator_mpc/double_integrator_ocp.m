function uopt = double_integrator_ocp(xt,N,P,Q,R,terminal_constraint)
if nargin < 6
    terminal_constraint = false;
end

cvx_begin
    variable u(N,1)
    variable x(2,N)
    f = cost(x,u,xt,P,Q,R,N);
    minimize f
    subject to
        for k = 0:N-1
            u(k+1) >= -0.5
            u(k+1) <= 0.5
%             if k < N-1
            x(:,k+1) >= [-5;-5]
            x(:,k+1) <= [5;5]
            if k == 0
                x(:,k+1) == double_integrator_dynamics(xt,u(k+1));
            else
                x(:,k+1) == double_integrator_dynamics(x(:,k),u(k+1))
            end
%             end
        end
cvx_end

if cvx_status == "Infeasible"
    uopt = [];
else
    uopt = u;
end

end


function f = cost(x,u,xt,P,Q,R,N)
f = 0;
for k = 0:N-1
    uk = u(k+1);
    if k==0
        xk = xt;
    else
        xk = x(:,k+1);
    end
    f = f + xk'*Q*xk + uk'*R*uk;
end

f = f + x(:,end)'*P*x(:,end);
end