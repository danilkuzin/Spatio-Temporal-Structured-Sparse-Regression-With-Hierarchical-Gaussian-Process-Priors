function [ F, precision, recall ] = F_measure_omega_spike( omega_est, omega_true)

T = size(omega_est, 2);
    
% vectorize
omega_true = omega_true(:);
omega_est = omega_est(:);
    

tp = sum((omega_true == 0) & (omega_est == 0));
fp = sum((omega_true == 1) & (omega_est == 0));
fn = sum((omega_true == 0) & (omega_est == 1));

precision = tp/(tp+fp);
recall = tp/(tp+fn);

if(precision+recall > 0)
    F = 2*(precision*recall/(precision+recall));
else
    F = 0;
end;

