function [converged, delta] = is_converged_mean_norm(Q, Q_old, epseps)

T = size(Q, 2);

beta = zeros(size(Q(1).beta.mean, 1), T);

beta_old = beta;

for t = 1 : T 
    beta_old(:, t) = Q_old(t).beta.mean;
    beta(:, t) = Q(t).beta.mean;
end

diff_beta = norm(beta(:) - beta_old(:), 'inf') / norm(beta_old(:), 'inf');

converged = diff_beta < epseps; %&& diff_gamma < epseps && diff_mu < epseps;

delta = max([diff_beta]);%, diff_gamma, diff_mu]);

end