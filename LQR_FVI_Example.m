% Initialization
format long
num_samples = 3000;
num_iterations = 10000;
w = ones(6, 1);  % Initial weights
threshold = 1e-3;
tau_max = 4;
Ts = 0.01;
discounted_f = 1;

% Cost functions
one_stage_cost = @(theta, omega, tau) 0.1 * theta^2 + 0.1 * omega^2 + tau^2;
J_approx = @(theta, omega, w) w' * [theta^2; theta * omega; omega^2; theta; omega; 1];

% Main loop
for iter = 1:num_iterations
    % writematrix(w, 'LQR_fixatone.dat', 'WriteMode', 'append');
    
    X = [];  % Feature matrix
    Y = [];  % Target cost vector
    
    for i = 1:num_samples
        % Sample state
        theta = 2 * pi * rand() - pi;
        omega = 4 * rand() - 2;
        
        % Compute features and target
        features = [theta^2; theta * omega; omega^2; theta; omega; 1];
        cost_fn = @(tau) one_stage_cost(theta, omega, tau) + discounted_f * J_approx(theta + Ts * omega, omega + Ts * tau, w);
        tau_opt = fminbnd(cost_fn, -tau_max, tau_max);
        target = cost_fn(tau_opt);
        
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
