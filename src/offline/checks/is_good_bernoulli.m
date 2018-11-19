function [res, msg] = is_good_bernoulli(factor)

MAX_NORMINVCDF_SUC_PROB = 1e+12;
assert(isfield(factor,'success_probability'));

res = true;
msg = '';

if any(isnan(factor.success_probability))
    res = false;
    msg = 'success probablity contains NaN values';
    return;
end

if any(abs(factor.success_probability) > MAX_NORMINVCDF_SUC_PROB)
    res = false;
    msg = 'norm cdfinv success probability values too large';
    return;
end