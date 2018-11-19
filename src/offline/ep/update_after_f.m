function [Q_updated] = update_after_f(Q, Qf, Qg, jitter, K, sigma2, X, T, run_checks)

Q_updated = Q;

for t=1:T
    Q_updated(t) = update_after_f_single_timestamp(Q(t), Qf(t), Qg(t), jitter, K, sigma2, X, run_checks);
end



