 function [results]=ep_func(Sigma, W, N, T, X, y, ...
     run_checks, sigma2_beta, K, sigma2, epseps, eta, eta_multiplier, ...
     available_omega, available_beta, eta_alls, max_iter_num, beta)

%consts
init_infty_inv = 1e-4;
init_infty = 1 / init_infty_inv;
jitter = 1e-9;
min_val = 1e-10;
max1_val = 1-(1e-10);
min_var = 1e-12;
max_var = 1e12;
infty_beta = 1e-2;
infty_gamma = 1e-6;


% EP
Qf_t = struct; % factor
Qf_t.beta = struct; % Gaussian
Qf_t.omega = struct; % Bernoulli
Qf_t.beta.natur_mean = zeros(N, 1); % Natural mean for Gaussian
Qf_t.beta.natur_var = init_infty_inv*ones(N, 1); % Natural variance for Gaussian, only diagonal elements
Qf_t.omega.success_probability = zeros(N, 1); % Success probability for Bernoulli

Qf = repmat(Qf_t, 1, T);

Qg_t = struct; % factor
Qg_t.beta = struct; % Gaussian
Qg_t.beta.natur_mean = zeros(N, 1); % Natural mean for Gaussian

Qg = repmat(Qg_t, 1, T);

for t = 1 : T
    Qg(t).beta.natur_mean = 1/sigma2 * X' * y(:, t);
end

Qh_t = struct; % factor
Qh_t.gamma = struct; % Gaussian
Qh_t.omega = struct; % Bernoulli
Qh_t.gamma.mean = zeros(N, 1); % Mean for Gaussian
Qh_t.omega.success_probability = zeros(N, 1); % Success probability for Bernoulli

Qh = repmat(Qh_t, 1, T);
Qh_joint = struct;
Qh_joint.gamma.common_var = init_infty*ones(N, 1);

Qr_t = struct; % factor
Qr_t.gamma = struct; % Gaussian
Qr_t.gamma.mean = zeros(N, 1);
Qr_t.mu = struct;
Qr_t.mu.mean = zeros(N, 1);

Qr = repmat(Qr_t, 1, T);
Qr_joint = struct;
Qr_joint.gamma.common_var = init_infty * ones(N, 1); % only diagonal elements are required
Qr_joint.mu.common_var = diag(init_infty * ones(N, 1));


Qu_t = struct;
Qu_t.mu_previous = struct;
Qu_t.mu_previous.mean = zeros(N, 1);
Qu_t.mu = struct;
Qu_t.mu.mean = zeros(N, 1);


Qu = repmat(Qu_t, 1, T);
Qu_joint = struct;
Qu_joint.mu.common_var = diag(init_infty * ones(N, 1));

Qu(1).mu_previous.mean = NaN;
Qu(1).mu.mean = NaN;


Q_t= struct; % joint distribution
Q_t.beta = struct; % Gaussian
Q_t.omega = struct; % Bernoulli
Q_t.beta.var = zeros(N, 1); % Variance for Gaussian
Q_t.beta.mean = zeros(N, 1); % Mean for Gaussian
Q_t.omega.success_probability = zeros(N, 1); % Success probability for Bernoulli

Q = repmat(Q_t, 1, T);
Q = update_after_f(Q, Qf, Qg, jitter, K, sigma2, X, T, run_checks);
Q_old = Q;


converged = false;
iter_num = 0;

delta_max = zeros(1, max_iter_num);

nmse = zeros(1, max_iter_num);
f_measure = zeros(1, max_iter_num);

while((~converged) && (iter_num < max_iter_num))
    iter_num = iter_num + 1;
    % f:
    Qf = update_f(Q, Qf, Qh, sigma2_beta, min_val, max1_val, min_var, max_var, infty_beta, eta, run_checks);
    Q = update_after_f(Q, Qf, Qg, jitter, K, sigma2, X, T, run_checks);
    
    % h:
    [Qh, Qh_joint] = update_h(Qf, Qh, Qh_joint, Qr, Qr_joint, min_var, max_var, infty_gamma, eta, run_checks);
   
    %r:
    [Qr, Qr_joint] = update_r(Qh, Qh_joint, Qr, Qr_joint, Qu, Qu_joint, Sigma, min_var, max_var, eta, run_checks);
    
    %u:
    [Qu, Qu_joint] = update_u(Qu, Qu_joint, Qr, Qr_joint, W, min_var, max_var, eta, jitter, run_checks);
    
    restored_beta = zeros(size(Q(1).beta.mean, 1), T);
    for t = 1 : T
        restored_beta(:, t) = Q(t).beta.mean;
    end    
    
    if (available_beta)
        nmse(iter_num) = compute_nmse(Q, beta);
    else 
        nmse(iter_num) = NMSE(X*restored_beta, y);
    end
    if (available_omega)
        Q = update_bern_param(Q, Qf, Qh);
        f_measure(iter_num) = compute_f_measure(Q, omega);
    end
    
    [converged, delta] = is_converged_mean_norm(Q,Q_old,epseps);
    delta_max(iter_num) = delta;
    
    eta_alls(iter_num) = eta;
    if (iter_num <= 2)
        converged = false;
    end

     if(converged)
         break;
     end
     
     Q_old = Q;
     eta = eta * eta_multiplier;    
end

Q = update_bern_param(Q, Qf, Qh);

results = struct;
restored_beta = zeros(N, T);
for t = 1 : T
    restored_beta(:, t) = Q(t).beta.mean;
end
results.restored_beta = restored_beta;
results.Q = Q;
results.Qf = Qf;
results.Qh = Qh;
results.Qh_joint = Qh_joint;
results.Qr = Qr;
results.Qr_joint = Qr_joint;
results.Qu = Qu;
results.Qu_joint = Qu_joint;
results.nmse = nmse;
results.iter_num = iter_num;
results.delta = delta_max;


