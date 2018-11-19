function [Qh_t_new, Qh_joint_new] = ...
        update_h_single_timestamp(Qf_t, Qh_t, Qh_joint, Qr_t, Qr_joint, min_var, max_var, infty_gamma, eta, T, run_checks)

Qh_t_new = Qh_t;
Qh_joint_new = Qh_joint;

if (run_checks)
    [res, msg] = is_good_gaussian(Qh_t.gamma);
    if (res == false)
        disp(msg);
        pause;
    end

    [res, msg] = is_good_bernoulli(Qf_t.omega);
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

gamma = struct;
gamma.cavity = struct;
gamma.cavity.var =  Qr_joint.gamma.common_var;
gamma.cavity.mean = Qr_t.gamma.mean;

omega = struct;
omega.cavity = struct;
omega.cavity.success_probability = Qf_t.omega.success_probability;

indexes_for_update = gamma.cavity.var > 0;
gamma.cavity.var(gamma.cavity.var <= 0) = 1;

a = gamma.cavity.mean ./ sqrt(1 + gamma.cavity.var); % temporary variable
Qh_new.omega.success_probability = a;

[log_phi_a, d_log_phi_a] = logphi(a);
log_phi_m_a = logphi(-a);
log_phi_p = logphi(omega.cavity.success_probability);
log_phi_m_p = logphi(-omega.cavity.success_probability);

log_z = logsum(log_phi_p + log_phi_a, log_phi_m_p + log_phi_m_a);

log_k1 = log(2) + log_phi_p + log_phi_a - log_z;
k1 = bound_values(exp(log_k1), min_var, max_var);

log_k2 = log_phi_a - log_z;
k2 = bound_values(exp(log_k2), min_var, max_var);

gamma.moments.first = gamma.cavity.mean + (k1 - k2) .* gamma.cavity.var .* d_log_phi_a ./ sqrt(1 + gamma.cavity.var);

log_k3 = log_phi_a + log_phi_p - log_z;
k3 = bound_values(exp(log_k3), min_var, max_var);

log_k4 = log_phi_m_p - log_z;
k4 = bound_values(exp(log_k4), min_var, max_var);

h = gamma.cavity.mean + d_log_phi_a .* gamma.cavity.var ./ sqrt(1 + gamma.cavity.var);
w = 2 * gamma.cavity.mean .* h - gamma.cavity.mean .^ 2 + gamma.cavity.var - d_log_phi_a .* a .* gamma.cavity.var .^ 2 ./ (1 + gamma.cavity.var);

gamma.moments.second = k4 .* (gamma.cavity.var + gamma.cavity.mean .^ 2 - normcdf(a) .* w) + k3 .* w;

% this is indeed Q capital new
Q_new.gamma.mean = gamma.moments.first;
Q_new.gamma.var = gamma.moments.second - gamma.moments.first .^ 2;
Q_new.gamma.var = bound_values(Q_new.gamma.var, min_var, max_var);

Qh_new.gamma.natur_var = 1 ./ Q_new.gamma.var - 1 ./ gamma.cavity.var;
Qh_new.gamma.natur_mean = Q_new.gamma.mean ./ Q_new.gamma.var - gamma.cavity.mean ./ gamma.cavity.var;

indexes_where_variance_is_negative = Qh_new.gamma.natur_var <= 0;
Qh_new.gamma.natur_var(indexes_where_variance_is_negative) = infty_gamma;
Q_new.gamma.var(indexes_where_variance_is_negative) = 1 ./ (Qh_new.gamma.natur_var(indexes_where_variance_is_negative) + 1 ./ gamma.cavity.var(indexes_where_variance_is_negative));
Qh_new.gamma.natur_mean(indexes_where_variance_is_negative) = Q_new.gamma.mean(indexes_where_variance_is_negative) ./ Q_new.gamma.var(indexes_where_variance_is_negative) - ...
    gamma.cavity.mean(indexes_where_variance_is_negative) ./ gamma.cavity.var(indexes_where_variance_is_negative);

%damping
Qh_joint_new.gamma.common_var(indexes_for_update) = 1 ./ ((eta / T) * Qh_new.gamma.natur_var(indexes_for_update) + ...
    (1 - eta / T) ./ Qh_joint.gamma.common_var(indexes_for_update));
Qh_joint_new.gamma.common_var = bound_values(Qh_joint_new.gamma.common_var, min_var, max_var);
Qh_t_new.gamma.mean(indexes_for_update) = Qh_joint_new.gamma.common_var(indexes_for_update) .* (eta * Qh_new.gamma.natur_mean(indexes_for_update) + ...
    (1 - eta) * Qh_t.gamma.mean(indexes_for_update) ./ Qh_joint.gamma.common_var(indexes_for_update));

Qh_t_new.omega.success_probability(indexes_for_update) = eta * Qh_new.omega.success_probability(indexes_for_update) + (1 - eta) * Qh_t.omega.success_probability(indexes_for_update);


if (run_checks)
    [res, msg] = is_good_gaussian(Qh_t_new.gamma);
    if (res == false)
        disp(msg);
        pause;
    end

    [res, msg] = is_good_gaussian(Qh_joint.gamma);
    if (res == false)
        disp(msg);
        pause;
    end
    
    
    [res, msg] = is_good_bernoulli(Qh_t_new.omega);
    if (res == false)
        disp(msg);
        pause;
    end
end
