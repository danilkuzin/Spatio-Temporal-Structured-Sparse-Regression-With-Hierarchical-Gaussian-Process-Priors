function nmse = compute_nmse(Q, beta)

T = size(Q, 2);

restored_beta = zeros(size(beta));
for t = 1 : T
    restored_beta(:, t) = Q(t).beta.mean;
end

nmse = NMSE(restored_beta, beta);

end