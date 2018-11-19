function [res, msg] = is_good_nonnatur_gaussian(factor)

MAX_VAR = 1e+12;
MAX_MEAN = 1e+12;

assert(isfield(factor,'var'));
assert(isfield(factor,'mean'));

res = true;
msg='';

if (isvector(factor.var))
    if any(factor.var <= 0)
        res = false;
        msg = 'var is not positive semidefinite';
        return;
    end
else
    [~,p] = chol(factor.var);
    if p>0
        res = false;
        msg = 'var is not positive semidefinite';
        return;
    end
end

if any(isnan(factor.var))
    res = false;
    msg = 'var contains NaN values';
    return;
end

if any(abs(factor.var) > MAX_VAR)
    res = false;
    msg = 'var contains large values';
    return;
end

if (~isvector(factor.var))
    if (~isequal(factor.var,factor.var'))
        res = false;
        msg = 'var is not symmetric';
        return;
    end
end

if any(isnan(factor.mean))
    res = false;
    msg = 'mean contains NaN values';
    return;
end

if any(abs(factor.mean) > MAX_MEAN)
    res = false;
    msg = 'mean contains large values';
    return;
end