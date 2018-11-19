function [Qr_t_new, Qr_joint_new] = update_r_single_timestamp_with_u1(t, Qh_t, Qh_joint, Qr_t, Qr_joint, Qu_t, Qu_t_plus_1, Qu_joint, Sigma, ...
    min_var, max_var, eta, T, run_checks)

Qr_t_new = Qr_t;
Qr_joint_new = Qr_joint;

if (run_checks)
    [res, msg] = is_good_gaussian(Qr_t.mu);
    if (res == false)
        disp(msg);
        pause;
    end
    
    [res, msg] = is_good_gaussian(Qr_t.gamma);
    if (res == false)
        disp(msg);
        pause;
    end

    [res, msg] = is_good_gaussian(Qh_t.gamma);
    if (res == false)
        disp(msg);
        pause;
    end
    
    [res, msg] = is_good_gaussian(Qr_joint.mu);
    if (res == false)
        disp(msg);
        pause;
    end
    
    [res, msg] = is_good_gaussian(Qr_joint.gamma);
    if (res == false)
        disp(msg);
        pause;
    end
    
    [res, msg] = is_good_gaussian(Qh_joint.gamma);
    if (res == false)
        disp(msg);
        pause;
    end
    
end

mu = struct;
mu.cavity = struct;

if t < T
    mu.cavity.var = 0.5 * diag(Qu_joint.mu.common_var);
    mu.cavity.mean = 0.5 * (Qu_t.mu.mean + Qu_t_plus_1.mu_previous.mean);
else
    mu.cavity.var = diag(Qu_joint.mu.common_var);
    mu.cavity.mean = Qu_t.mu.mean;
end

gamma = struct;
gamma.cavity = struct;
gamma.cavity.var = Qh_joint.gamma.common_var;
gamma.cavity.mean = Qh_t.gamma.mean;


Qr_new.gamma.var = mu.cavity.var + diag(Sigma);
Qr_new.gamma.mean = mu.cavity.mean;

Qr_new.mu.var = diag(gamma.cavity.var) + Sigma;
Qr_new.mu.mean = gamma.cavity.mean;


% damping
Qr_joint_new.gamma.common_var = 1 ./ ((eta / T) ./ Qr_new.gamma.var + (1 - eta / T) ./ Qr_joint.gamma.common_var);
Qr_joint_new.gamma.common_var = bound_values(Qr_joint_new.gamma.common_var, min_var, max_var);
Qr_t_new.gamma.mean = Qr_joint_new.gamma.common_var .* (eta * Qr_new.gamma.mean ./ Qr_new.gamma.var + (1 - eta) * Qr_t.gamma.mean ./ Qr_joint.gamma.common_var);

Qr_t_new.mu.mean = eta * Qr_new.mu.mean + (1 - eta) * Qr_t.mu.mean;
Qr_joint_new.mu.common_var = (eta / T) * Qr_new.mu.var + (1 - eta / T) * Qr_joint.mu.common_var;
Qr_joint_new.mu.common_var = bound_values_matrix(Qr_joint_new.mu.common_var, min_var, max_var);


if (run_checks)
    [res, msg] = is_good_gaussian(Qr_t_new.gamma);
    if (res == false)
        disp(msg);
        pause;
    end
    
    [res, msg] = is_good_gaussian(Qr_joint.gamma);
    if (res == false)
        disp(msg);
        pause;
    end

    [res, msg] = is_good_gaussian(Qr_t_new.mu);
    if (res == false)
        disp(msg);
        pause;
    end
    
    [res, msg] = is_good_gaussian(Qr_joint.mu);
    if (res == false)
        disp(msg);
        pause;
    end
end