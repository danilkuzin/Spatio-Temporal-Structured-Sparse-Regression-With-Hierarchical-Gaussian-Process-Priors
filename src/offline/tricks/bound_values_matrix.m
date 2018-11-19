function [ output ] = bound_values_matrix( matrix, lower, upper, run_checks )

assert(size(matrix,1)==size(matrix,2));
output = matrix;

indexes_to_bound = find(diag(output) <= lower);
for i=1:numel(indexes_to_bound)
    ind = indexes_to_bound(i);
    output(ind,:)=0;
    output(:,ind)=0;
    output(ind,ind)=lower;
    if run_checks
        disp('bounded matrix value lower');
    end
end

indexes_to_bound = find(diag(output) >= upper);
for i=1:numel(indexes_to_bound)
    ind = indexes_to_bound(i);
    output(ind,ind)=upper;
    if run_checks
        disp('bounded matrix value upper');
    end
end