function Q_updated_t = update_after_f_single_timestamp(Q_t, Qf_t, Qg_t, jitter, K, sigma2, X, run_checks)

Q_updated_t = Q_t;

if (run_checks)
    [res, msg] = is_good_gaussian(Qf_t.beta);
    if (res == false)
        disp(msg);
        pause;
    end

%     [res, msg] = is_good_gaussian(Qg_t.beta);
%     if (res == false)
%         disp(msg);
%         pause;
%     end
end

I = (1 + jitter) * eye(K);
theta_scaled = sigma2 * Qf_t.beta.natur_var.^2;
beta = struct;
beta.natur_mean = Qg_t.beta.natur_mean + Qf_t.beta.natur_mean;
B = bsxfun(@rdivide, X, sqrt(Qf_t.beta.natur_var)');
L = chol((B*B.')/sigma2 + I, 'lower');
R = L \ X;
R_sqr = sum(R.^2,1)';
Q_updated_t.beta.var = 1./Qf_t.beta.natur_var - R_sqr./theta_scaled;
beta.natur_mean_scaled = beta.natur_mean ./ Qf_t.beta.natur_var;
Q_updated_t.beta.mean = beta.natur_mean_scaled - R.'*(R*beta.natur_mean_scaled)./(sigma2.*Qf_t.beta.natur_var);

if (run_checks)
    [res, msg] = is_good_nonnatur_gaussian(Q_updated_t.beta);
    if (res == false)
        disp(msg);
        pause;
    end
end