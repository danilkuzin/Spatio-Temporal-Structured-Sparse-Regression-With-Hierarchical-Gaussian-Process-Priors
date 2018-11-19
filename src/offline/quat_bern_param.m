function result = quat_bern_param(divident, divider)

result = norminv((normcdf(divident) - normcdf(divident) * normcdf(divider))/(normcdf(divident)+normcdf(divider)-2*normcdf(divident)*normcdf(divider)), 0 , 1);

end