function Qf_t_new = update_f_single_timestamp(Q_t, Qf_t, Qh_t, sigma2_beta, ...
    min_val, max1_val, min_var, max_var, infty_beta, eta, run_checks)

Qf_t_new = Qf_t;
Qf_new = Qf_t;

if (run_checks)
    [res, msg] = is_good_nonnatur_gaussian(Q_t.beta);
    if (res == false)
        disp(msg);
        pause;
    end

    [res, msg] = is_good_gaussian(Qf_t.beta);
    if (res == false)
        disp(msg);
        pause;
    end

    [res, msg] = is_good_bernoulli(Qh_t.omega);
    if (res == false)
        disp(msg);
        pause;
    end
end

beta = struct;
beta.cavity = struct;
beta.cavity.natur_mean = Q_t.beta.mean ./ Q_t.beta.var - Qf_t.beta.natur_mean;
beta.cavity.natur_var = 1 ./ Q_t.beta.var - Qf_t.beta.natur_var;
indexes_for_update = beta.cavity.natur_var > 0;
beta.cavity.natur_var(beta.cavity.natur_var <= 0) = 1;

beta.cavity.mean = beta.cavity.natur_mean ./ beta.cavity.natur_var;
beta.cavity.var = 1 ./ beta.cavity.natur_var;

omega = struct;
omega.cavity = struct;
omega.cavity.success_probability = Qh_t.omega.success_probability;

r2 = lognormpdf(0, beta.cavity.mean, beta.cavity.var + sigma2_beta) - lognormpdf(0, beta.cavity.mean, beta.cavity.var); % temporary variable
q2 = 1 ./ (1 + exp(r2)); % temporary variable
Qf_new.omega.success_probability = norminv(min(max(q2,min_val),max1_val));

omega.cavity.logqout = logphi(omega.cavity.success_probability) - logphi(-omega.cavity.success_probability);
C = bound_values(1 + exp(omega.cavity.logqout - r2), min_var, max_var); % temporary variable
D = beta.cavity.mean .* beta.cavity.natur_var .* sigma2_beta ./ (1 + sigma2_beta .* beta.cavity.natur_var);  % temporary variable
beta.moments.first = D ./ C;
beta.moments.second = (D .^ 2 + sigma2_beta ./ (1 + sigma2_beta .* beta.cavity.natur_var)) ./ C;
beta.mean = beta.moments.first;
beta.var = beta.moments.second - beta.moments.first .^ 2;
Qf_new.beta.natur_var = 1 ./ beta.var - beta.cavity.natur_var;
Qf_new.beta.natur_mean = beta.mean ./ beta.var - beta.cavity.natur_mean;

indexes_where_variance_is_negative = Qf_new.beta.natur_var <= 0;
% assign var to positive infinity
Qf_new.beta.natur_var(indexes_where_variance_is_negative) = infty_beta;
Qf_new.beta.natur_var = bound_values(Qf_new.beta.natur_var, min_var, max_var);
beta.var(indexes_where_variance_is_negative) = 1 ./ (Qf_new.beta.natur_var(indexes_where_variance_is_negative) + beta.cavity.natur_var(indexes_where_variance_is_negative));
Qf_new.beta.natur_mean(indexes_where_variance_is_negative) = beta.mean(indexes_where_variance_is_negative) ./ beta.var(indexes_where_variance_is_negative) - beta.cavity.natur_mean(indexes_where_variance_is_negative);

%damping
Qf_t_new.beta.natur_var(indexes_for_update) = eta * Qf_new.beta.natur_var(indexes_for_update) + (1 - eta) * Qf_t.beta.natur_var(indexes_for_update);
Qf_t_new.beta.natur_mean(indexes_for_update) = eta * Qf_new.beta.natur_mean(indexes_for_update) + (1 - eta) * Qf_t.beta.natur_mean(indexes_for_update);
%TODO consider replacing with prod_bern param
Qf_t_new.omega.success_probability(indexes_for_update) = eta * Qf_new.omega.success_probability(indexes_for_update) + (1 - eta) * Qf_t.omega.success_probability(indexes_for_update);

if (run_checks)
    [res, msg] = is_good_gaussian(Qf_t_new.beta);
    if (res == false)
        disp(msg);
        pause;
    end

    [res, msg] = is_good_bernoulli(Qf_t_new.omega);
    if (res == false)
        disp(msg);
        pause;
    end
end