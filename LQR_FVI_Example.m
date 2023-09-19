% Initialization
format long
num_samples = 3000;
num_iterations = 10000;
w = ones(6, 1);  % Initial weights
threshold = 1e-3;
u_max = 4;
Ts = 0.01;
discounted_f = 1;

% Cost functions
one_stage_cost = @(x0, x1, u) 0.1 * x0^2 + 0.1 * x1^2 + u^2;
J_approx = @(x0, x1, w) w' * [x0^2; x0 * x1; x1^2; x0; x1; 1];

% Main loop
for iter = 1:num_iterations
    % writematrix(w, 'LQR_fixatone.dat', 'WriteMode', 'append');
    
    X = [];  % Feature matrix
    Y = [];  % Target cost vector
    
    for i = 1:num_samples
        % Sample state
        x0 = 4 * rand() - 2;
        x1 = 4 * rand() - 2;
        
        % Compute features and target
        features = [x0^2; x0 * x1; x1^2; x0; x1; 1];
        cost_fn = @(u) one_stage_cost(x0, x1, u) + discounted_f * J_approx(x0 + Ts * x1, x1 + Ts * u, w);
        u_opt = fminbnd(cost_fn, -u_max, u_max);
        target = cost_fn(u_opt);
        
        % Store data
        X = [X; features'];
        Y = [Y; target];
    end
    
    % Update weights
    w_new = X \ Y;
    
    % Convergence check
    if norm(w_new - w) < threshold
        disp(iter);
        break;
    end
    
    w = w_new;
    disp(w);
end

% Final output
S = [w(2), 0.5 * w(1); 0.5 * w(1), w(3)];
disp(S);
disp(w);
