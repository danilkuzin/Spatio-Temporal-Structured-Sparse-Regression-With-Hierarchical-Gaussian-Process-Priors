% script runs offline ep inference for the sparse spatio-temporal
% structured model on the sythetic data (to reproduce the results from the
% paper)
diary on;

root_path = '../';
src_path = [root_path, 'src/offline/'];
addpath(genpath(src_path));

SData = load([root_path, 'data/synthetic/', num2str(1), '/data_',num2str(5),'_size.mat']);

[N, T] = size(SData.beta);

sigma2_beta = 10000;
sigma2 = 2;
epseps = 1e-3;
eta = 1-1e-3;
eta_multiplier = 1-1e-4;
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
        W(i,j) = a1 * exp(-(i-j)^2/(2*s1^2));
    end
end

W = (W + W')/2;
W = W + jitter * eye(size(W));


Sigma = zeros(N, N);
for i=1:N
    for j=1:N
        Sigma(i,j) = a2 * exp(-(i-j)^2/(2*s2^2));
    end
end

Sigma = (Sigma + Sigma') / 2;
Sigma = Sigma + jitter * eye(size(Sigma));


max_iter_num = 1000;

available_omega = false;
available_beta = true;

run_checks = false;

for experiment_number=1:10
    for size_of_y=10:5:55
        disp(['experiment: ', num2str(experiment_number), ' y size: ', num2str(size_of_y)]);
        data_path = [root_path, 'data/synthetic/', num2str(experiment_number), '/data_',num2str(size_of_y),'_size.mat'];
        results_path = [root_path, 'results/offline/synthetic/', ...
            num2str(experiment_number), '/', num2str(size_of_y), '/'];

        SData = load(data_path);

        K = size(SData.y, 1);
        eta_alls = zeros(1, max_iter_num);

        tic
        results = ep_func(Sigma, W, N, T, SData.X, SData.y, ...
             run_checks, sigma2_beta, K, sigma2, epseps, eta, eta_multiplier, ...
             available_omega, available_beta, eta_alls, max_iter_num, SData.beta);
        results.time = toc; 
        
        mkdir(results_path)
        save([results_path, 'results.mat'], 'results', '-v7.3');
    end
end
rmpath(src_path);

disp('finish!')
