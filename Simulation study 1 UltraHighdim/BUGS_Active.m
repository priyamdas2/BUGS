clearvars;
addpath('./BUGS/');
addpath('./Data/');
addpath('./supp funs/');

% Scenario settings
scenario = 1;
num_reps = 10;
n = 200;
p = 100000;
% Choose replicate

%% BUGS options

opts.niter = 5000;
opts.burnin = 1000;
opts.thin = 5;
opts.display_each = 100;
opts.verbose = true;

opts.tau0 = 0.001;           % smaller => stronger global shrinkage
opts.sigma_eta = 0.002;      % guidance useful => put higher
opts.Cz = 2;
opts.a_c = 10;   % E(c^2) =b_c/(a_c-1), small c^2 => more shrinkage, small a_c  more heavy tail
opts.b_c = 100;
opts.a_sigma = 2;
opts.b_sigma = 2;
use_beta_summary = 'mean';
%use_beta_summary = 'median';
sel_thresh = 0.01;

%% BUGS-Active specific options
opts.save_full_beta = true;      % needed for out2.beta_mean / out2.beta_median
opts.max_p_full_save = p;        % must be >= p
opts.save_top_beta = 20000;        % 100 for low dims

% active-set controls
opts.active_beta_thresh = 0.001;
opts.active_cap = 25000;            % 1000 for low dims
opts.min_active = min(p,20000);      % 100 for low dims
opts.guidance_keep = min(p,20000);   % 200 for low dims
opts.lambda2_inactive = 1e-6;



parfor rep = 1:num_reps    
    
    fprintf('performing data rep: %d | scenario: %d | n: %d | p: %d\n', rep, scenario, n, p);
    
    % File names
    fileX = sprintf('Data_Scenario_%d_X_n_%d_p_%d_rep_%d.csv',scenario, n, p, rep);
    fileY = sprintf('Data_Scenario_%d_Y_n_%d_p_%d_rep_%d.csv',scenario, n, p, rep);
    fileB = sprintf('Data_Scenario_%d_beta_n_%d_p_%d.csv',scenario, n, p);
    
    % Read data
    X = readmatrix(fileX);
    y_raw = readmatrix(fileY);
    beta_true = readmatrix(fileB);
    
    % Force correct shapes
    y_raw = y_raw(:);
    beta_true = beta_true(:);
    
    % Basic checks
    if size(X,1) ~= n || size(X,2) ~= p
        error('X has dimension %d x %d, expected %d x %d.', size(X,1), size(X,2), n, p);
    end
    if length(y_raw) ~= n
        error('y has length %d, expected %d.', length(y_raw), n);
    end
    if length(beta_true) ~= p
        error('beta_true has length %d, expected %d.', length(beta_true), p);
    end
    
    % Standardize y and corresponding beta_true
    y_raw = y_raw - mean(y_raw);
    sy = std(y_raw);
    
    if sy == 0
        error('Standard deviation of y_raw is zero.');
    end
    
    y = y_raw / sy;
    beta_true_std = beta_true / sy;
    
    
    
    
    % ----------------------------
    % Fit BUGS and record time
    % ----------------------------
    tic;
    out = BUGS_Active_mcmc(X, y, opts);
    comp_time = toc;
    
    % ----------------------------
    % Choose beta estimate
    % ----------------------------
    switch lower(use_beta_summary)
        case 'mean'
            beta_hat = out.beta_mean;
        case 'median'
            beta_hat = out.beta_median;
        otherwise
            error('use_beta_summary must be ''mean'' or ''median''.');
    end
    
    % ----------------------------
    % Compute metrics
    % ----------------------------
    metrics = compute_metrics_from_beta_BUGS(X, beta_true_std, beta_hat, comp_time, sel_thresh);
    
    % ----------------------------
    % Save one-row results for this rep
    % ----------------------------
    
    
    T = struct2table(metrics);
    
    
    outfile_rep = sprintf('Results/Output_scen_%d_n_%d_p_%d_BUGSactive_rep_%d.csv', scenario, n, p, rep);
    writetable(T, outfile_rep);
    
    disp(T);
    fprintf('Saved metrics to: %s\n', outfile_rep);
    
end