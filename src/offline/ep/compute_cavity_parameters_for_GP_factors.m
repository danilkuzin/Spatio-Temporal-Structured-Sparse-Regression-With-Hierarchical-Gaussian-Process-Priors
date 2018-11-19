function [cavity_var, cavity_mean] = compute_cavity_parameters_for_GP_factors(var_first, var_second, mean_first, mean_second, jitter)
    N = numel(mean_first);
    L = chol(var_first + var_second + jitter * eye(N), 'lower');
    R = L \ var_first;
    cavity_var = var_first - R' * R; 
    
    cavity_mean = mean_first - (var_first / L') * (L \ mean_first) + ...
        mean_second - (var_second / L') * (L \ mean_second);
end