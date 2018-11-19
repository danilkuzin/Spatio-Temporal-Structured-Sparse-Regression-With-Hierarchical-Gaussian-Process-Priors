T = 50;
N = 100;
sparsity = 0.95;
num_of_groups = 2;
start_column = 1;
sqrt_sigma2_beta = 100;
min_non_zero_elements = 2;

number_of_experiments = 10;

for experiment_number = 1:number_of_experiments
    
    rng(1000 + experiment_number, 'twister');
    
    omega = ones(N, T);
    slabs = randsample(N,num_of_groups);
    slabs = sort(slabs);
    
    group_sizes = poissrnd(0.5 * (1-sparsity)*N/num_of_groups, num_of_groups, 1);
    
    start_group = zeros(num_of_groups, 1);
    end_group = zeros(num_of_groups, 1);
    
    group_disappeared = false(num_of_groups, 1);
    
    for i=1:num_of_groups
        start_group(i) = slabs(i);
        end_group(i) = min(slabs(i) + group_sizes(i), N);
        omega(start_group(i):end_group(i), start_column) = 0;
    end
    
    for t=start_column+1:T
        evolutions = datasample([1:3], 2*num_of_groups, 'Weights', [0.5*sparsity, 1, 0.5*(1 / sparsity)]);
        for i=1:num_of_groups
            if (evolutions(2*i-1)==1)
                start_group(i) = max(start_group(i)-1, 1);
            elseif (evolutions(2*i-1)==2)
                start_group(i) = start_group(i);
            elseif (evolutions(2*i-1)==3)
                start_group(i) = min(start_group(i)+1, N);
            end
            
            if (evolutions(2*i)==3)
                end_group(i) = max(end_group(i)-1, 1);
            elseif (evolutions(2*i)==2)
                end_group(i) = end_group(i);
            elseif (evolutions(2*i)==1)
                end_group(i) = min(end_group(i)+1, N);
            end
            
            if end_group(i) > start_group(i)
                if end_group(i) - start_group(i) > 0.5 * (1 - sparsity) * N
                    end_group(i) = max(end_group(i)-1, 1);
                end
            end
            
            if (start_group(i) + min_non_zero_elements > end_group(i))
                if all(group_disappeared(setdiff(1:num_of_groups, i)) == true)
                    end_group(i) = min(start_group(i)+min_non_zero_elements, N);
                end
            end
            if start_group(i) > end_group(i)
                group_disappeared(i) = true;
            else
                group_disappeared(i) = false;
            end    
            
            omega(start_group(i):end_group(i), t) = 0;
        end
    end
        
    beta = zeros(N, T);
    
    for t = 1: T
        for i = 1: N
            if omega(i, t) == 0
                beta(i, t) = randn(1) * sqrt_sigma2_beta;
            end
        end
    end
    
    figure;
    surf(beta);
    
    y_sizes = 5:5:100;
    
    mkdir(['../../data/synthetic/', num2str(experiment_number)]);
    for j=y_sizes
        X = randn(j, N);
        y = X*beta;
        save(['../../data/synthetic/', num2str(experiment_number), '/data_', num2str(j), '_size.mat'], 'y', 'X', 'beta', 'omega');
    end
    
end
