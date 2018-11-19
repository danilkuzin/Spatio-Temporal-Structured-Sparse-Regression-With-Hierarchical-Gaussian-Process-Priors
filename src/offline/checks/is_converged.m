function [converged, delta_max] = is_converged(Q, Q_old, epseps, T)
    converged = false;
    delta_max = 0;
    for t = 1 : T
        delta = max([compute_change(Q_old(t).beta.mean, Q(t).beta.mean), compute_change(diag(Q_old(t).beta.var), diag(Q(t).beta.var)), ...
                compute_change(Q_old(t).gamma.natur_mean, Q(t).gamma.natur_mean), ...
                compute_change(diag(Q_old(t).gamma.natur_var), diag(Q(t).gamma.natur_var)), ...
                compute_change(Q_old(t).omega.success_probability, Q(t).omega.success_probability), ...
                compute_change(Q_old(t).mu.natur_mean, Q(t).mu.natur_mean), ...
                compute_change(diag(Q_old(t).mu.natur_var), diag(Q(t).mu.natur_var))]);
        if max(delta) > delta_max
            delta_max = max(delta);
        end
    end
    
    if (delta_max < epseps)
        converged = true;
    end
end