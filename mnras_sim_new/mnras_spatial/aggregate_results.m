function [x, residual, asnr, ssnr, asnr_log, ssnr_log, acpu, scpu, arun, srun, total_cpu_time, total_runtime, iteration_number] = ...
    aggregate_results(algo, filename, ncores_data, ncores_prior, x0, squared_operator_norm, Q, upsilon0)
%%
% Produce the images and metrics reported in the MNRAS paper
% ``A Faceted Prior for Scalable Wideband Imaging: Application to Radio
% Astronomy''
%
%-------------------------------------------------------------------------%
%%
% Author: P.-A. Thouvenin.
% Last modified: [../../2020]
%-------------------------------------------------------------------------%
%%
% NOTES:
%
%
% TODO: load operator norm, compute snr (print in a .txt file), compute
% timing, relate to number of CPUs
%
%-------------------------------------------------------------------------%
%%
N = [size(x0, 1), size(x0, 2)];
nChannels = size(x0, 3);
x = zeros(N(1), N(2), nChannels);
residual = zeros(N(1), N(2), nChannels);

% Reconstruction metrics
% DR = @(x, res, n) sqrt(prod(N))*squeeze(max(max(x,[],2),[],1))*n./norm2D(res); % per band input
norm2D = @(x) squeeze(sqrt(sum(sum(x.^2, 2), 1)));
SNR = @(x, x0) 20*log10(norm2D(x0)./norm2D(x - x0));
% SNR_log = @(x, x0) 20*log10(norm2D(log10(x0+eps))./norm2D(log10(x+eps)-log10(x0+eps)));
SNR_log = @(x, x0) 20*log10(norm2D(log10(1 + x0./reshape(upsilon0, [1, 1, numel(upsilon0)])))./norm2D(log10(1 + x./reshape(upsilon0, [1, 1, numel(upsilon0)]))-log10(1 + x0./reshape(upsilon0, [1, 1, numel(upsilon0)]))));

% asnr = zeros(numel(Qx), 1);
% asnr_log = zeros(numel(Qx), 1);
% vsnr = zeros(numel(Qx), 1);
% vsnr_log = zeros(numel(Qx), 1);
% iteration_number = zeros(numel(Qx), 1);
% total_runtime = zeros(numel(Qx), 1);  % total runtime (in s)
% total_cpu_time = zeros(numel(Qx), 1); % total CPU time (in s)
% aruntime = zeros(numel(Qx), 1);       % average runtime (per iteration, in s)
% vruntime = zeros(numel(Qx), 1);       % variance runtime
% acpu_time = zeros(numel(Qx), 1);      % average cpu time (per iter., in s)
% vcpu_time = zeros(numel(Qx), 1);      % variance cpu time
% atime_facet = zeros(numel(Qx), 1);    % average time (facet, per iter., in s)
% atime_data = zeros(numel(Qx), 1);     % variance
% vtime_facet = zeros(numel(Qx), 1);    % average time (data, per iter., in s)
% vtime_data = zeros(numel(Qx), 1);     % variance

switch algo
    case 'sara'
    
        total_runtime = 0;  % assuming all sub-cubes are running in parallel
        total_cpu_time = 0; %
        arun = 0;

        atime_l11 = 0; % average time (per iteration)
        atime_data = 0;
        atime_master = 0;
        vtime_l11 = 0; % variance
        vtime_data = 0;
        vtime_master = 0;
        iteration_number = 0;

        sum_runtime_sqr = 0;
        sum_l11_sqr = 0;
        sum_master_sqr = 0;
        sum_data_sqr = 0;
        sum_cpu_sqr = 0;

        for l = 1:nChannels
            load(filename(l), 'xsol', 't_l11', 't_master', 't_data', 'end_iter', 'res')
            x(:,:,l) = xsol;
            residual(:,:,l) = res/squared_operator_norm(l);
            
            % compute mean and variance over all the files? just 1 for now
            arun = arun + sum(end_iter(end_iter > 0)); % average runtime per iteration, over all sub-problems
            sum_runtime_sqr = sum_runtime_sqr + sum(end_iter(end_iter > 0).^2);
            
            atime_l11 = atime_l11 + sum(t_l11(t_l11 > 0));
            atime_master = atime_master + sum(t_master(t_master > 0));
            atime_data = atime_data + sum(t_data(t_data > 0));
            sum_l11_sqr = sum_l11_sqr + (ncores_prior^2)*sum(t_l11(t_l11 > 0).^2);
            sum_master_sqr = sum_master_sqr + sum(t_master(t_master > 0).^2);
            sum_data_sqr = sum_data_sqr + sum(t_data(t_data > 0).^2);
            
            % average number of iterations over all sub-problems
            iteration_number = iteration_number + sum(end_iter > 0);
            total_runtime = max(total_runtime, sum(end_iter(end_iter > 0)));
            total_cpu_time = total_cpu_time + ncores_prior*sum(t_l11(t_l11 > 0)) + ...
                + sum(t_master(t_master > 0)) ...
                + sum(t_data(t_data > 0));
            sum_cpu_sqr = sum_cpu_sqr + sum((ncores_prior*t_l11(t_l11 > 0) + t_master(t_master > 0) ...
                + t_data(t_data > 0)).^2);
        end
        
        arun = arun/iteration_number;
        % atime_l11 = atime_l11/iteration_number;
        % atime_master = atime_master/iteration_number;
        % atime_data = atime_data/iteration_number;
        %
        srun = sqrt((sum_runtime_sqr - iteration_number*arun^2)/(iteration_number - 1));
        % stime_l11 = sqrt((sum_l11_sqr - iteration_number*atime_l11^2)/(iteration_number - 1));
        % stime_master = sqrt((sum_master_sqr - iteration_number*atime_master^2)/(iteration_number - 1));
        % stime_data = sqrt((sum_data_sqr - iteration_number*atime_data^2)/(iteration_number - 1));
        %
        acpu = total_cpu_time/iteration_number;
        scpu = sqrt((sum_cpu_sqr - iteration_number*acpu^2)/(iteration_number - 1));
        iteration_number = round(iteration_number/nChannels); % average number of iterations for each SARA problem
        
        % compute SNR
        a = SNR(x, x0);
        asnr = mean(a);
        ssnr = std(a);
        a = SNR_log(x, x0);
        asnr_log = mean(a);
        ssnr_log = std(a);

        ncpu = nChannels*(ncores_data + ncores_prior + 1);

    case 'cw'

        load(filename, 'xsol', 'res', 't_facet', 't_data', 'end_iter')
        x = xsol;
        residual = res./squared_operator_norm;
        clear xsol res;

        % compute SNR
        a = SNR(x, x0);
        asnr = mean(a);
        ssnr = std(a);
        a = SNR_log(x, x0);
        asnr_log = mean(a);
        ssnr_log = std(a);

        % timing
        iteration_number = sum((end_iter > 0));
        %
        total_runtime = sum(end_iter(end_iter > 0));
        arun = total_runtime/iteration_number; % average runtime per iteration
        srun = std(end_iter(end_iter > 0));
        %
        atime_facet = mean(t_facet(t_facet > 0));          % t_facet already averaged over the Q facets
        stime_facet = sqrt((Q^4)*var(t_facet(t_facet > 0))); % multiplicative factor to account for the initial averaging
        atime_data = mean(t_data(t_data > 0));
        stime_data = sqrt((ncores_data^2)*var(t_data(t_data > 0)));
        %
        a_ = Q^2*t_facet(t_facet > 0) + ncores_data*t_data(t_data > 0);
        acpu = Q^2*atime_facet + ncores_data*atime_data;
        scpu = std(a_);
        total_cpu_time = sum(a_);

        ncpu = min(ncores_prior + ncores_data + 1, 36);

    case 'hypersara'
        load(filename, 'xsol', 'res', 't_l21', 't_nuclear', 't_master', 't_data', 'end_iter')
        x = xsol;
        residual = res./squared_operator_norm;
        clear xsol res;

        % compute SNR
        a = SNR(x, x0);
        asnr = mean(a);
        ssnr = std(a);
        a = SNR_log(x, x0);
        asnr_log = mean(a);
        ssnr_log = std(a);

        % timing
        a_ = t_l21(t_l21 > 0) + t_nuclear(t_nuclear > 0) + t_master(t_master > 0);
        b_ = t_data(t_data > 0);
        iteration_number = sum((end_iter > 0));
        %
        atime_facet = mean(a_);
        stime_facet = std(a_);
        atime_data = mean(b_);
        stime_data = std(b_);
        %
        total_runtime = sum(end_iter(end_iter > 0));
        arun = total_runtime/iteration_number; % average runtime per iteration
        srun = std(end_iter(end_iter > 0));
        %
        a_ = a_ + ncores_data*b_;
        total_cpu_time = sum(a_);
        acpu = atime_facet + ncores_data*atime_data; % average cpu time per iteration
        scpu = std(a_);

        ncpu = min(ncores_prior + ncores_data + 1, 36);
end

%% Display results (table)
fprintf("%s, asnr = %2.2f (%1.2e), asnr_log = %2.2f (%1.2e), cpu_cores = %i, iteration_number = %i \n", ...
    algo, asnr, ssnr, asnr_log, ssnr_log, ncpu, iteration_number)
fprintf(" arun (s) = %.2f (%1.2e), total_runtime (h) = %2.2f \n", ...   
    arun, srun, total_runtime/3600);
fprintf("acpu (s) = %.2f (%1.2e) , total_cpu_time (h) = %2.2f \n", ...
    acpu, scpu, total_cpu_time/3600);

end