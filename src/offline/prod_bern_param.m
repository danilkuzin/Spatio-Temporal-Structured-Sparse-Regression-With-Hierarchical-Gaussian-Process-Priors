function result = prod_bern_param(mult_first, mult_second)

l = numel(mult_first);
result = zeros(size(mult_first));

for i = 1 : l
    result(i) = norminv(normcdf(mult_first(i))*normcdf(mult_second(i))/(1-normcdf(mult_first(i))-normcdf(mult_second(i))+2*normcdf(mult_first(i))*normcdf(mult_second(i))), 0 , 1);
end

end