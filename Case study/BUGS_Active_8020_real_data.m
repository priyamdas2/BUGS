% Takes ~ 50 sec to load data
clearvars -except X_raw_FULL y_raw_FULL
clc;
addpath('./BUGS/');
addpath('./Dummy data/');
addpath('./Age methylation data/');

% output folder
outdir = './Output/';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Keep only first keep_n_upto rows and first keep_p_upto columns
% Use 'all' to keep everything
keep_n_upto = 'all';   % e.g. 200, 500, 'all'    / max 1051
keep_p_upto = 'all';   % e.g. 1000, 50000, 'all'  / max 865859

p = 865859;

%% =========================
%  BUGS options
%  =========================
opts.niter = 1000;
opts.burnin = 200;
opts.thin = 5;
opts.display_every = 2;
opts.verbose = true;

opts.tau0 = 0.001;           
opts.sigma_eta = 1e-6;      
opts.Cz = 2;
opts.a_c = 10;              
opts.b_c = 100;
opts.a_sigma = 2;
opts.b_sigma = 2;

%% =========================
%  BUGS-Active specific options
%  =========================
opts.save_full_beta = true;      % needed for beta_mean / beta_median
opts.max_p_full_save = p;        % must be >= p
opts.save_top_beta = min(p, 200);          % for low-dimensional dummy data, save all top entries

% active-set controls
opts.active_beta_thresh = 0.001;
opts.active_cap = min(p, 500);
opts.min_active = min(p, 200);
opts.guidance_keep = min(p, 200);
opts.lambda2_inactive = 1e-3;



%% =========================
%  User settings
%  =========================
seed = 1;
rng(seed);

train_frac = 0.80;
use_beta_summary = 'mean';   % 'mean' or 'median'
%% =========================
%  File names
%  =========================
fileX = 'X_matrix_final.mat';
fileY = 'Pheno_for_MATLAB.csv';

%% =========================
%  Read full data (only if not already loaded)
%  =========================

% ---- X ----
if ~exist('X_raw_FULL','var') || isempty(X_raw_FULL)
    fprintf('Loading full X...\n');
    t_load = tic;
    S = load(fileX);
    fn = fieldnames(S);
    X_raw_FULL = S.(fn{1});
    clear S fn;
    fprintf('Time to load full X: %.2f seconds\n', toc(t_load));
else
    fprintf('Full X already loaded. Skipping load.\n');
end

% ---- Y ----
if ~exist('y_raw_FULL','var') || isempty(y_raw_FULL)
    Ty = readtable(fileY);
    y_raw_FULL = Ty.Age_months;
    y_raw_FULL = y_raw_FULL(:);
    clear Ty;
else
    fprintf('Full y already loaded. Skipping load.\n');
end

%% =========================
%  Infer full dimensions
%  =========================
[n_full, p_full] = size(X_raw_FULL);

fprintf('Loaded full data with n = %d and p = %d.\n', n_full, p_full);

%% =========================
%  Basic checks on full data
%  =========================
if length(y_raw_FULL) ~= n_full
    error('y has length %d, but X has %d rows.', length(y_raw_FULL), n_full);
end

%% =========================
%  Apply keep_n_upto / keep_p_upto
%  =========================
if ischar(keep_n_upto) || isstring(keep_n_upto)
    if strcmpi(keep_n_upto, 'all')
        n_keep = n_full;
    else
        error('keep_n_upto must be numeric or ''all''.');
    end
else
    n_keep = min(keep_n_upto, n_full);
end

if ischar(keep_p_upto) || isstring(keep_p_upto)
    if strcmpi(keep_p_upto, 'all')
        p_keep = p_full;
    else
        error('keep_p_upto must be numeric or ''all''.');
    end
else
    p_keep = min(keep_p_upto, p_full);
end

X_raw = X_raw_FULL(1:n_keep, 1:p_keep);
y_raw = y_raw_FULL(1:n_keep);

[n, p] = size(X_raw);

fprintf('Using working data with n = %d and p = %d.\n', n, p);

%% =========================
%  80-20 split
%  =========================
perm = randperm(n);
n_train = floor(train_frac * n);
idx_train = perm(1:n_train);
idx_test  = perm(n_train+1:end);

X_train_raw = X_raw(idx_train, :);
X_test_raw  = X_raw(idx_test, :);

y_train_raw = y_raw(idx_train);
y_test_raw  = y_raw(idx_test);

%% =========================
%  Standardize X using training stats only
%  =========================
muX = mean(X_train_raw, 1);
sdX = std(X_train_raw, 0, 1);

% protect against zero-variance columns
sdX(sdX == 0) = 1;

X_train = (X_train_raw - muX) ./ sdX;
X_test  = (X_test_raw  - muX) ./ sdX;

%% =========================
%  Standardize y using training stats only
%  =========================
muY = mean(y_train_raw);
sdY = std(y_train_raw);

if sdY == 0
    error('Standard deviation of training y is zero.');
end

y_train = (y_train_raw - muY) / sdY;
y_test  = (y_test_raw  - muY) / sdY;

%% =========================
%  Converting to double for optimized BLAS
%  =========================

X_train = double(X_train);
X_test  = double(X_test);
y_train = double(y_train);
y_test  = double(y_test);

%% =========================
%  Fit BUGS_Active_mcmc on training data
%  =========================
tic;
out = BUGS_Active_mcmc(X_train, y_train, opts);
comp_time = toc;

%% =========================
%  Choose beta estimate
%  =========================
switch lower(use_beta_summary)
    case 'mean'
        beta_hat = out.beta_mean;
    case 'median'
        beta_hat = out.beta_median;
    otherwise
        error('use_beta_summary must be ''mean'' or ''median''.');
end

beta_hat = beta_hat(:);

if length(beta_hat) ~= p
    error('beta_hat has length %d, expected %d.', length(beta_hat), p);
end

%% =========================
%  Test prediction
%  =========================
yhat_test_std = X_test * beta_hat;
resid_test_std = y_test - yhat_test_std;

% back-transform to original y scale
yhat_test_raw = muY + sdY * yhat_test_std;
resid_test_raw = y_test_raw - yhat_test_raw;

%% =========================
%  Prediction metrics
%  =========================
rmse_test = sqrt(mean((y_test_raw - yhat_test_raw).^2));
mae_test  = mean(abs(y_test_raw - yhat_test_raw));

if std(y_test_raw) > 0 && std(yhat_test_raw) > 0
    corr_test = corr(y_test_raw, yhat_test_raw);
else
    corr_test = NaN;
end

sst = sum((y_test_raw - mean(y_test_raw)).^2);
sse = sum((y_test_raw - yhat_test_raw).^2);
if sst > 0
    r2_test = 1 - sse / sst;
else
    r2_test = NaN;
end

%% =========================
%  Selected CpGs count (simple threshold on chosen beta summary)
%  =========================
sel_thresh = 0.01;
n_selected = sum(abs(beta_hat) > sel_thresh);

%% =========================
%  Save test predictions for predicted-vs-observed plot
%  =========================
Tpred = table();
Tpred.obs_index = idx_test(:);
Tpred.y_true_raw = y_test_raw(:);
Tpred.y_pred_raw = yhat_test_raw(:);
Tpred.resid_raw = resid_test_raw(:);
Tpred.y_true_std = y_test(:);
Tpred.y_pred_std = yhat_test_std(:);
Tpred.resid_std = resid_test_std(:);

outfile_pred = fullfile(outdir, 'Output_8020_guided_cv_test_predictions.csv');
writetable(Tpred, outfile_pred);

%% =========================
%  Save one-row metrics table
%  =========================
Tmetrics = table();
Tmetrics.n = n;
Tmetrics.p = p;
Tmetrics.n_train = n_train;
Tmetrics.n_test = n - n_train;
Tmetrics.beta_summary = string(use_beta_summary);
Tmetrics.rmse_test = rmse_test;
Tmetrics.mae_test = mae_test;
Tmetrics.corr_test = corr_test;
Tmetrics.r2_test = r2_test;
Tmetrics.n_selected = n_selected;
Tmetrics.comp_time_sec = comp_time;
Tmetrics.y_train_mean = muY;
Tmetrics.y_train_sd = sdY;

outfile_metrics = fullfile(outdir, 'Output_8020_guided_cv_metrics.csv');
writetable(Tmetrics, outfile_metrics);

%% =========================
%  Save hyperparameter draws
%  =========================
Thyper_draws = table();

if isfield(out, 'sigma2_draws')
    Thyper_draws.sigma2 = out.sigma2_draws(:);
end
if isfield(out, 'tau_draws')
    Thyper_draws.tau = out.tau_draws(:);
end
if isfield(out, 'c2_draws')
    Thyper_draws.c2 = out.c2_draws(:);
end
if isfield(out, 'eta_draws')
    Thyper_draws.eta = out.eta_draws(:);
end
if isfield(out, 'active_size_draws')
    Thyper_draws.active_size = out.active_size_draws(:);
end

outfile_hdraws = fullfile(outdir, 'Output_8020_guided_cv_hyper_draws.csv');
writetable(Thyper_draws, outfile_hdraws);

%% =========================
%  Save hyperparameter summaries
%  =========================
param_names = {};
post_mean = [];
post_median = [];
post_sd = [];
ci_2p5 = [];
ci_97p5 = [];

if isfield(out, 'sigma2_draws')
    x = out.sigma2_draws(:);
    param_names{end+1,1} = 'sigma2';
    post_mean(end+1,1) = mean(x);
    post_median(end+1,1) = median(x);
    post_sd(end+1,1) = std(x);
    ci_2p5(end+1,1) = prctile(x, 2.5);
    ci_97p5(end+1,1) = prctile(x, 97.5);
end

if isfield(out, 'tau_draws')
    x = out.tau_draws(:);
    param_names{end+1,1} = 'tau';
    post_mean(end+1,1) = mean(x);
    post_median(end+1,1) = median(x);
    post_sd(end+1,1) = std(x);
    ci_2p5(end+1,1) = prctile(x, 2.5);
    ci_97p5(end+1,1) = prctile(x, 97.5);
end

if isfield(out, 'c2_draws')
    x = out.c2_draws(:);
    param_names{end+1,1} = 'c2';
    post_mean(end+1,1) = mean(x);
    post_median(end+1,1) = median(x);
    post_sd(end+1,1) = std(x);
    ci_2p5(end+1,1) = prctile(x, 2.5);
    ci_97p5(end+1,1) = prctile(x, 97.5);
end

if isfield(out, 'eta_draws')
    x = out.eta_draws(:);
    param_names{end+1,1} = 'eta';
    post_mean(end+1,1) = mean(x);
    post_median(end+1,1) = median(x);
    post_sd(end+1,1) = std(x);
    ci_2p5(end+1,1) = prctile(x, 2.5);
    ci_97p5(end+1,1) = prctile(x, 97.5);
end

if isfield(out, 'active_size_draws')
    x = out.active_size_draws(:);
    param_names{end+1,1} = 'active_size';
    post_mean(end+1,1) = mean(x);
    post_median(end+1,1) = median(x);
    post_sd(end+1,1) = std(x);
    ci_2p5(end+1,1) = prctile(x, 2.5);
    ci_97p5(end+1,1) = prctile(x, 97.5);
end

Thyper_summary = table(param_names, post_mean, post_median, post_sd, ci_2p5, ci_97p5, ...
    'VariableNames', {'parameter','mean','median','sd','ci_2p5','ci_97p5'});

outfile_hsummary = fullfile(outdir, 'Output_8020_guided_cv_hyper_summary.csv');
writetable(Thyper_summary, outfile_hsummary);

%% =========================
%  Save train/test split info
%  =========================
Tsplit = table();
Tsplit.obs_index = (1:n)';
Tsplit.in_train = ismember((1:n)', idx_train);
Tsplit.in_test  = ismember((1:n)', idx_test);

outfile_split = fullfile(outdir, 'Output_8020_guided_cv_split.csv');
writetable(Tsplit, outfile_split);

%% =========================
%  Console output
%  =========================
disp('====================================================');
disp('BUGS_Active_mcmc dummy 80-20 run completed.');
disp(Tmetrics);
disp(Thyper_summary);
fprintf('Saved files to folder: %s\n', outdir);
disp('====================================================');