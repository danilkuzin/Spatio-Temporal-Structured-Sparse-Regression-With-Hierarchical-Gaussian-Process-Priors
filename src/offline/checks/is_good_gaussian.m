function [res, msg] = is_good_gaussian(factor)

MAX_VAR = 1e+12;
MAX_MEAN = 1e+12;

res = true;
msg='';

if isfield(factor, 'natur_var')
    if (isvector(factor.natur_var))
        if any(factor.natur_var <= 0)
            res = false;
            msg = 'natur var is not positive semidefinite';
            return;
        end
    else
        [~,p] = chol(factor.natur_var);
        if p>0
            res = false;
            msg = 'natur var is not positive semidefinite';
            return;
        end
    end

    if any(isnan(factor.natur_var))
        res = false;
        msg = 'natur var contains NaN values';
        return;
    end

    if (~isvector(factor.natur_var))
        if (~isequal(factor.natur_var,factor.natur_var'))
            res = false;
            msg = 'natur var is not symmetric';
            return;
        end
    end

    if any(abs(factor.natur_var) > MAX_VAR)
        res = false;
        msg = 'natur var contains large values';
        return;
    end
end

if isfield(factor, 'common_natur_var')
    if (isvector(factor.common_natur_var))
        if any(factor.common_natur_var <= 0)
            res = false;
            msg = 'natur var is not positive semidefinite';
            return;
        end
    else
        [~,p] = chol(factor.common_natur_var);
        if p>0
            res = false;
            msg = 'natur var is not positive semidefinite';
            return;
        end
    end

    if any(isnan(factor.common_natur_var))
        res = false;
        msg = 'natur var contains NaN values';
        return;
    end

    if (~isvector(factor.common_natur_var))
        if (~isequal(factor.natur_var,factor.natur_var'))
            res = false;
            msg = 'natur var is not symmetric';
            return;
        end
    end

    if any(abs(factor.common_natur_var) > MAX_VAR)
        res = false;
        msg = 'natur var contains large values';
        return;
    end
end


if isfield(factor,'natur_mean')
    if any(isnan(factor.natur_mean))
        res = false;
        msg = 'natur mean contains NaN values';
        return;
    end

    if any(abs(factor.natur_mean) > MAX_MEAN)
        res = false;
        msg = 'natur mean contains large values';
        return;
    end
end