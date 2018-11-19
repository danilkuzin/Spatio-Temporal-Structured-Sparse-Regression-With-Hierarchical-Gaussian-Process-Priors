function results = run_online(beta, X, y, eta, eta_multiplier)
% function that performs online ep inference for the sparse 
% spatio-temporal structured model on given data: beta, X, and y. It 
% specifies hyperparameters, generates covariance matrices and run 
% inference.     
    [N, T] = size(beta);
    K = size(y, 1);
    
    sigma2_beta = 10000;
    sigma2 = 100;
    epseps = 1e-3;
    scalar_for_mu0 = 10;
    mu0 = scalar_for_mu0 * ones(N, 1);
    s1 = 15;
    s2 = 10;
    a1 = 10;
    a2 = 10;
    
    jitter = 1e-9;
    
    W=zeros(N,N);
    for i=1:N
        for j=1:N
            W(i,j) = a1 * exp(-(i-j)^2/(2*s1^2)); % 1 * ...
        end
    end
    
    W = (W + W')/2;
    W = W + jitter * eye(size(W));
    
    Sigma = zeros(N, N);
    for i=1:N
        for j=1:N
            Sigma(i,j) = a2 * exp(-(i-j)^2/(2*s2^2)); % 0.2 * ...
        end
    end
    
    Sigma = (Sigma + Sigma') / 2;
    Sigma = Sigma + jitter * eye(size(Sigma));
    
    
    max_iter_num = 1000;
    
    available_omega = false;
    available_beta = true;
    eta_alls = zeros(1, max_iter_num);
    
    run_checks = false;
    
    [results] = offline_and_online_ep_inference(Sigma, W, N, T, X, y, ...
        run_checks, sigma2_beta, K, sigma2, epseps, eta, eta_multiplier, mu0, ...
        available_omega, available_beta, eta_alls, max_iter_num, beta);
    
end