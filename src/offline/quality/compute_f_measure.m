function f_measure = compute_f_measure(Q, omega)

threshold = 0.5;

T = size(Q, 2);

restored_omega_success_probability = zeros(size(omega));
for t = 1 : T
    restored_omega_success_probability(:, t) = normcdf(Q(t).omega.success_probability);
end

restored_omega = restored_omega_success_probability > threshold;

f_measure = F_measure_omega_spike(restored_omega, omega);

end