clearvars;

addpath('./BUGS/');
addpath('./Data/');
addpath('./supp funs/');

% ============================================================
% Scenario settings
% ============================================================

scenario = 2;
num_reps = 10;

n = 200;
p = 10000;

% ============================================================
% SIS-RHS options
% ============================================================

opts.niter = 5000;
opts.burnin = 1000;
opts.thin = 5;

% Regularized horseshoe hyperparameters
% Keep these identical to the corresponding BUGS settings
% for a fair comparison.
opts.tau0 = 0.001;

opts.a_c = 10;
opts.b_c = 100;

opts.a_sigma = 2;
opts.b_sigma = 2;

% Slice sampler settings
opts.slice_w = 0.5;
opts.slice_m = 100;

% Output / display settings
opts.verbose = true;
opts.display_every = 100;
opts.save_lambda_max = 100;

% Posterior beta summary
use_beta_summary = 'mean';
% use_beta_summary = 'median';

% Threshold used for variable-selection metrics
sel_thresh = 0.01;

% Create Results folder if needed
if ~exist('Results', 'dir')
    mkdir('Results');
end


% ============================================================
% Loop over simulation replicates
% ============================================================

parfor rep = 1:num_reps
    
    rng(rep);
    
    fprintf(['performing data rep: %d | scenario: %d | ' ...
             'n: %d | p: %d\n'], ...
             rep, scenario, n, p);
    
    
    % ========================================================
    % File names
    % ========================================================
    
    fileX = sprintf( ...
        'Data_Scenario_%d_X_n_%d_p_%d_rep_%d.csv', ...
        scenario, n, p, rep);
    
    fileY = sprintf( ...
        'Data_Scenario_%d_Y_n_%d_p_%d_rep_%d.csv', ...
        scenario, n, p, rep);
    
    fileB = sprintf( ...
        'Data_Scenario_%d_beta_n_%d_p_%d.csv', ...
        scenario, n, p);
    
    
    % ========================================================
    % Read data
    % ========================================================
    
    X = readmatrix(fileX);
    y_raw = readmatrix(fileY);
    beta_true = readmatrix(fileB);
    
    % Force correct shapes
    y_raw = y_raw(:);
    beta_true = beta_true(:);
    
    
    % ========================================================
    % Basic checks
    % ========================================================
    
    if size(X,1) ~= n || size(X,2) ~= p
        error(['X has dimension %d x %d, expected ' ...
               '%d x %d.'], ...
               size(X,1), size(X,2), n, p);
    end
    
    if length(y_raw) ~= n
        error('y has length %d, expected %d.', ...
              length(y_raw), n);
    end
    
    if length(beta_true) ~= p
        error('beta_true has length %d, expected %d.', ...
              length(beta_true), p);
    end
    
    
    % ========================================================
    % Standardize y and corresponding beta_true
    % ========================================================
    
    y_raw = y_raw - mean(y_raw);
    sy = std(y_raw);
    
    if sy == 0
        error('Standard deviation of y_raw is zero.');
    end
    
    y = y_raw / sy;
    beta_true_std = beta_true / sy;
    
    
    % ========================================================
    % SIS + Regularized Horseshoe
    %
    % Computational time includes BOTH:
    %   1. SIS screening
    %   2. RHS MCMC
    % ========================================================
    
    tic;
    
    
    % --------------------------------------------------------
    % Step 1: Sure Independence Screening
    %
    % Fan and Lv (2008):
    % retain d = floor(n / log(n)) variables with largest
    % absolute marginal correlations with y.
    % --------------------------------------------------------
    
    d_sis = min(p, floor(n / log(n)));
    
    % Predictor standard deviations
    x_sd = std(X, 0, 1);
    
    % Since y is centered, X' * y is proportional to the
    % centered marginal covariance.
    %
    % Dividing by sd(X_j) gives the same ranking as
    % absolute marginal correlations because sd(y) is common
    % to every predictor.
    sis_score = abs(X' * y) ./ x_sd';
    
    % Handle zero-variance or invalid predictors
    sis_score(~isfinite(sis_score)) = -Inf;
    
    % Rank predictors by SIS score
    [~, sis_order] = sort(sis_score, 'descend');
    
    % Keep top d_sis predictors
    sis_idx = sis_order(1:d_sis);
    
    % Reduced design matrix
    X_sis = X(:, sis_idx);
    
    
    % --------------------------------------------------------
    % Step 2: Regularized Horseshoe on SIS-selected variables
    % --------------------------------------------------------
    
    % IMPORTANT:
    % BUGSRHS_mcmc resets the RNG internally using opts.seed.
    % Therefore give each replicate its own seed.
    opts_rep = opts;
    opts_rep.seed = rep;
    
    out = BUGSRHS_mcmc(X_sis, y, opts_rep);
    
    
    % Total SIS + RHS computation time
    comp_time = toc;
    
    
    % ========================================================
    % Choose posterior beta estimate in reduced space
    % ========================================================
    
    switch lower(use_beta_summary)
        
        case 'mean'
            beta_hat_reduced = out.beta_mean;
            
        case 'median'
            beta_hat_reduced = out.beta_median;
            
        otherwise
            error(['use_beta_summary must be ' ...
                   '''mean'' or ''median''.']);
    end
    
    beta_hat_reduced = beta_hat_reduced(:);
    
    
    % ========================================================
    % Expand estimate back to original p-dimensional space
    %
    % Predictors excluded by SIS are assigned coefficient 0.
    % Therefore TPR/FPR/FDR/MCC are calculated relative to
    % all original p predictors.
    % ========================================================
    
    beta_hat = zeros(p, 1);
    beta_hat(sis_idx) = beta_hat_reduced;
    
    
    % ========================================================
    % Compute metrics
    % ========================================================
    
    metrics = compute_metrics_from_beta_BUGS( ...
        X, ...
        beta_true_std, ...
        beta_hat, ...
        comp_time, ...
        sel_thresh);
    
    
    % ========================================================
    % Save one-row results for this replicate
    % ========================================================
    
    T = struct2table(metrics);
    
    outfile_rep = sprintf( ...
        ['Results/Output_scen_%d_n_%d_p_%d_' ...
         'SIS_RHS_rep_%d.csv'], ...
         scenario, n, p, rep);
    
    writetable(T, outfile_rep);
    
    
    % ========================================================
    % Display
    % ========================================================
    
    fprintf(['SIS retained %d of %d predictors ' ...
             'for replicate %d.\n'], ...
             d_sis, p, rep);
    
    disp(T);
    
    fprintf('Saved metrics to: %s\n', outfile_rep);
    
end