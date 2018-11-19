function [Qu_t_new, Qu_joint_new] = update_u_single_timestamp(t, Qu_t, Qu_tm1, Qu_t_plus_1, Qu_joint, ...
    Qr_t, Qr_tm1, Qr_joint, W, min_var, max_var, eta, T, jitter, run_checks)

Qu_t_new = Qu_t;
Qu_joint_new = Qu_joint;

if (run_checks)
    [res, msg] = is_good_gaussian(Qu_t.mu);
    if (res == false)
        disp(msg);
        pause;
    end
    
    [res, msg] = is_good_gaussian(Qu_joint.mu);
    if (res == false)
        disp(msg);
        pause;
    end
end
    
mu = struct;
mu.cavity = struct;

if t < T
    [mu.cavity.var, mu.cavity.mean] = ...
        compute_cavity_parameters_for_GP_factors(Qr_joint.mu.common_var, Qu_joint.mu.common_var, Qr_t.mu.mean, Qu_t_plus_1.mu_previous.mean, jitter);
else
    mu.cavity.var = Qr_joint.mu.common_var;
    mu.cavity.mean = Qr_t.mu.mean;
end


mu_previous = struct;
mu_previous.cavity = struct;

if t == 2
    mu_previous.cavity.var = Qr_joint.mu.common_var;
    mu_previous.cavity.mean = Qr_tm1.mu.mean;
else
    [mu_previous.cavity.var, mu_previous.cavity.mean] = ...
        compute_cavity_parameters_for_GP_factors(Qr_joint.mu.common_var, Qu_joint.mu.common_var, Qr_tm1.mu.mean, Qu_tm1.mu.mean, jitter);
end


Qu_new.mu.var = mu_previous.cavity.var + W;
Qu_new.mu.mean = mu_previous.cavity.mean;

Qu_new.mu_previous.var = mu.cavity.var + W;
Qu_new.mu_previous.mean = mu.cavity.mean;

% damping
Qu_t_new.mu_previous.mean = eta * Qu_new.mu_previous.mean + ...
    (1 - eta) * Qu_t.mu_previous.mean;
Qu_joint_new.mu.common_var = ...
    (eta / T / 2) * Qu_new.mu_previous.var + ...
    (1 - eta / T / 2) * Qu_joint.mu.common_var;
Qu_joint_new.mu.common_var = bound_values_matrix(Qu_joint_new.mu.common_var, min_var, max_var, run_checks);


Qu_t_new.mu.mean = eta * Qu_new.mu.mean + (1 - eta) * Qu_t.mu.mean;
Qu_joint_new.mu.common_var = (eta / T / 2) * Qu_new.mu.var + ...
    (1 - eta / T / 2) * Qu_joint_new.mu.common_var;
Qu_joint_new.mu.common_var = bound_values_matrix(Qu_joint_new.mu.common_var, min_var, max_var, run_checks);

if (run_checks)
    [res, msg] = is_good_gaussian(Qu_t_new.mu);
    if (res == false)
        disp(msg);
        pause;
    end

    [res, msg] = is_good_gaussian(Qu_t_new.mu_previous);
    if (res == false)
        disp(msg);
        pause;
    end
    
    [res, msg] = is_good_gaussian(Qu_joint.mu);
    if (res == false)
        disp(msg);
        pause;
    end
end