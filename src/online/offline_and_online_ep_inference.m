function [final_results] = offline_and_online_ep_inference(Sigma, W, N, T, X, y, ...
            run_checks, sigma2_beta, K, sigma2, epseps, eta, eta_multiplier, ...
            available_omega, available_beta, eta_alls, max_iter_num, beta, ...
            offline_training_size, batch_size)
    
    % parameters
    jitter = 1e-9;
    
    if (T <= offline_training_size)
        error('T size');
    end
    batches_amount = floor((T - offline_training_size)/batch_size);
    
    %prepare offline training data
    y_offline = y(:, 1:offline_training_size);
    if available_beta
        beta_offline = beta(:, 1:offline_training_size);
    else
        beta_offline = [];
    end
        
    %disp('running offline training')
    [offline_results] = ep_func(Sigma, W, N, offline_training_size, X, y_offline, ...
             run_checks, sigma2_beta, K, sigma2, epseps, eta, eta_multiplier, ...
             available_omega, available_beta, eta_alls, max_iter_num, beta_offline);  
    
    restored_beta = offline_results.restored_beta;
    
    [mu_posterior.var, mu_posterior.mean] = ...
        compute_cavity_parameters_for_GP_factors(offline_results.Qr_joint.mu.common_var, ...
                                                 offline_results.Qu_joint.mu.common_var, ...
                                                 offline_results.Qr(offline_training_size).mu.mean, ...
                                                 offline_results.Qu(offline_training_size).mu.mean, jitter);    
        
    for i=1:batches_amount
        beta_prior.var = sigma2_beta * ones(N, 1);
        
        %prepare batch training data
        y_batch = y(:, offline_training_size+(i-1)*batch_size + 1 : offline_training_size+i*batch_size);
        if available_beta
            beta_batch = beta(:, offline_training_size+(i-1)*batch_size + 1 : offline_training_size+i*batch_size);
        else
            beta_batch = [];
        end
    
        %predict        
        mu_prior_predict.var = mu_posterior.var + W;
        mu_prior_predict.mean = mu_posterior.mean;
    
        %update
        [online_results] = ep_online_func(Sigma, ...
            W, N, batch_size, X, y_batch, ...
            run_checks, beta_prior.var, K, sigma2, epseps, eta, eta_multiplier, ...
            mu_prior_predict.mean, mu_prior_predict.var, ...
            available_omega, available_beta, eta_alls, max_iter_num, beta_batch);
        
        restored_beta(:, ...
            offline_training_size+(i-1)*batch_size + 1 : ...
                offline_training_size+i*batch_size) = online_results.restored_beta;
            
        [mu_posterior.var, mu_posterior.mean] = ...
            compute_cavity_parameters_for_GP_factors(online_results.Qr_joint.mu.common_var, ...
                                                 online_results.Qu_joint.mu.common_var, ...
                                                 online_results.Qr(batch_size).mu.mean, ...
                                                 online_results.Qu(batch_size).mu.mean, jitter);
    end
    
    final_results.restored_beta = restored_beta;
end