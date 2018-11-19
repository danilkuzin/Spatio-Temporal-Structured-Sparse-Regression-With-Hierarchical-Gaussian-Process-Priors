function [Qu, Qu_joint] = update_u(Qu, Qu_joint, Qr, Qr_joint, W, min_var, max_var, eta, jitter, run_checks)

T = size(Qu, 2);

for t=2:T-1
    [Qu(t), Qu_joint] = update_u_single_timestamp(t, Qu(t), Qu(t-1), Qu(t+1), Qu_joint, ...
        Qr(t), Qr(t-1), Qr_joint, W, min_var, max_var, eta, T, jitter, run_checks); 
end


[Qu(T), Qu_joint] = update_u_single_timestamp(T, Qu(T), Qu(T-1), [], Qu_joint, ...
    Qr(T), Qr(T-1), Qr_joint, W, min_var, max_var, eta, T, jitter, run_checks); 


