clearvars;
addpath('./BUGS/');
addpath('./Data/');
addpath('./supp funs/');

%% =========================
%  Simulation settings
%  =========================
sensitivity_casenum = 2;% Sensitivity-analysis case selector: 1 to 12


scenario = 1;             
num_reps = 10;
n = 200;
p = 1000;


% Posterior summary used for point estimate
use_beta_summary = 'mean';

% Selection threshold for support recovery metrics
sel_thresh = 0.01;

%% =========================
%  MCMC controls
%  =========================
base_opts = struct();
base_opts.niter = 5000;
base_opts.burnin = 1000;
base_opts.thin = 5;
base_opts.display_every = 100;   
base_opts.verbose = true;

% Get hyperparameters for the selected sensitivity case
[opts, case_label] = get_sensitivity_opts(base_opts, sensitivity_casenum);

fprintf('Running sensitivity case %d: %s\n', sensitivity_casenum, case_label);
disp(opts);

%% =========================
%  Main loop
%  =========================
parfor rep = 1:num_reps
    rng(rep);

    fprintf(['performing data rep: %d | scenario: %d | n: %d | p: %d ' ...
             '| sensitivity case: %d\n'], ...
             rep, scenario, n, p, sensitivity_casenum);

    % ----------------------------
    % File names for data
    % ----------------------------
    fileX = sprintf('Data_Scenario_%d_X_n_%d_p_%d_rep_%d.csv', scenario, n, p, rep);
    fileY = sprintf('Data_Scenario_%d_Y_n_%d_p_%d_rep_%d.csv', scenario, n, p, rep);
    fileB = sprintf('Data_Scenario_%d_beta_n_%d_p_%d.csv', scenario, n, p);

    % ----------------------------
    % Read data
    % ----------------------------
    X = readmatrix(fileX);
    y_raw = readmatrix(fileY);
    beta_true = readmatrix(fileB);

    % Force correct shapes
    y_raw = y_raw(:);
    beta_true = beta_true(:);

    % ----------------------------
    % Basic checks
    % ----------------------------
    if size(X,1) ~= n || size(X,2) ~= p
        error('X has dimension %d x %d, expected %d x %d.', ...
            size(X,1), size(X,2), n, p);
    end

    if length(y_raw) ~= n
        error('y has length %d, expected %d.', length(y_raw), n);
    end

    if length(beta_true) ~= p
        error('beta_true has length %d, expected %d.', length(beta_true), p);
    end

    % ----------------------------
    % Standardize y and corresponding beta_true
    % ----------------------------
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
    out = BUGS_mcmc(X, y, opts);
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
    % Add case metadata to results
    % ----------------------------
    metrics.sensitivity_casenum = sensitivity_casenum;
    metrics.case_label = string(case_label);

    metrics.tau0 = opts.tau0;
    metrics.sigma_eta = opts.sigma_eta;
    metrics.Cz = opts.Cz;
    metrics.a_c = opts.a_c;
    metrics.b_c = opts.b_c;
    metrics.a_sigma = opts.a_sigma;
    metrics.b_sigma = opts.b_sigma;

    metrics.niter = opts.niter;
    metrics.burnin = opts.burnin;
    metrics.thin = opts.thin;

    % ----------------------------
    % Save one-row results for this rep
    % ----------------------------
    T = struct2table(metrics);

    outfile_rep = sprintf(['Results/Output_casenum_%d_n_%d_p_%d_' ...
                           'BUGS_rep_%d.csv'], ...
                           sensitivity_casenum, n, p, rep);

    writetable(T, outfile_rep);

    disp(T);
    fprintf('Saved metrics to: %s\n', outfile_rep);
end