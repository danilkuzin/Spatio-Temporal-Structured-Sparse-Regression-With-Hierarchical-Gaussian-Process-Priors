function Q_new = update_bern_param(Q, Qf, Qh)

T = size(Q, 2);

Q_new = Q;

for t = 1 : T
    Q_new(t).omega.success_probability = prod_bern_param(Qf(t).omega.success_probability, Qh(t).omega.success_probability);
end

end