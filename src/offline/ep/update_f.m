function Qf_new = update_f(Q, Qf, Qh, sigma2_beta, min_val, max1_val, min_var, max_var, infty_beta, eta, run_checks)

T = size(Q, 2);

Qf_new = Qf;

%TODO consider replacing with parfor, or (:) vector operations
for t=1:T
    Qf_new(t) = update_f_single_timestamp(Q(t), Qf(t), Qh(t), sigma2_beta, ...
        min_val, max1_val, min_var, max_var, infty_beta, eta, run_checks);
end