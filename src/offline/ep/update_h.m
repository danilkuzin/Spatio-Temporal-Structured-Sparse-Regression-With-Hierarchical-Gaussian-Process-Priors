function [Qh, Qh_joint] = update_h(Qf, Qh, Qh_joint, Qr, Qr_joint, min_var, max_var, infty_gamma, eta, run_checks)

T = size(Qh, 2);

for t=1:T
    [Qh(t), Qh_joint] = update_h_single_timestamp(Qf(t), Qh(t), Qh_joint, Qr(t), Qr_joint, min_var, max_var, infty_gamma, eta, T, run_checks);
end