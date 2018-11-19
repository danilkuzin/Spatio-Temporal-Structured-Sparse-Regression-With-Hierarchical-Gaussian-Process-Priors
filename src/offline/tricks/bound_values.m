function [ output ] = bound_values( values, lower, upper )
% The function takes an input array 'values' and bounds the values between
% 'lower' and 'upper' by thresholding.
%
% Danil Kuzin, Olga Isupova
    output = min(max(values, lower), upper);
end