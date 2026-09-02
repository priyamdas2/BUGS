% ============================================================
% REAL DATA 80-20 ANALYSIS:
%   SIS-LASSO
%   SIS-SCAD
%   SIS-RHS
%
% All three methods:
%   - use exactly the same 80/20 split
%   - use training-only standardization
%   - use exactly the same SIS-selected predictors
%   - save separate method-specific outputs
%
% Requires:
%   BUGSRHS_mcmc.m in ./BUGS/
%   Statistics and Machine Learning Toolbox for MATLAB lasso()
% ============================================================

clearvars -except X_raw_FULL y_raw_FULL
clc;

addpath('./BUGS/');
addpath('./Dummy data/');
addpath('./Age methylation data/');

% ============================================================
% Output folder
% ============================================================

outdir = './Output/';

if ~exist(outdir, 'dir')
    mkdir(outdir);
end


% ============================================================
% Keep rows / columns
% ============================================================

% Use 'all' to keep everything
keep_n_upto = 'all';   % max 1051
keep_p_upto = 'all';   % max 865859


% ============================================================
% User settings
% ============================================================

seed = 1;
rng(seed);

train_frac = 0.80;

% Selection threshold used consistently across methods
sel_thresh = 0.01;

% CV settings
nfolds = 10;

% SCAD settings
scad_gamma = 3.7;
scad_num_lambda = 100;
scad_lambda_ratio = 1e-3;
scad_max_iter = 1000;
scad_tol = 1e-6;

% RHS posterior summary
use_beta_summary = 'mean';
% use_beta_summary = 'median';


% ============================================================
% RHS MCMC options
% ============================================================

opts_rhs.niter = 1000;
opts_rhs.burnin = 200;
opts_rhs.thin = 5;

opts_rhs.tau0 = 0.001;

opts_rhs.a_c = 10;
opts_rhs.b_c = 100;

opts_rhs.a_sigma = 2;
opts_rhs.b_sigma = 2;

opts_rhs.slice_w = 0.5;
opts_rhs.slice_m = 100;

opts_rhs.verbose = true;
opts_rhs.display_every = 100;
opts_rhs.save_lambda_max = 100;

opts_rhs.seed = seed;


% ============================================================
% File names
% ============================================================

fileX = 'X_matrix_final.mat';
fileY = 'Pheno_for_MATLAB.csv';


% ============================================================
% Read full data only if not already loaded
% ============================================================

% ----------------------------
% X
% ----------------------------

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


% ----------------------------
% y
% ----------------------------

if ~exist('y_raw_FULL','var') || isempty(y_raw_FULL)

    Ty = readtable(fileY);

    y_raw_FULL = Ty.Age_months;
    y_raw_FULL = y_raw_FULL(:);

    clear Ty;

else

    fprintf('Full y already loaded. Skipping load.\n');

end


% ============================================================
% Full dimensions
% ============================================================

[n_full, p_full] = size(X_raw_FULL);

fprintf('Loaded full data with n = %d and p = %d.\n', ...
    n_full, p_full);


% ============================================================
% Basic check
% ============================================================

if length(y_raw_FULL) ~= n_full
    error('y has length %d, but X has %d rows.', ...
        length(y_raw_FULL), n_full);
end


% ============================================================
% Apply keep_n_upto / keep_p_upto
% ============================================================

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

fprintf('Using working data with n = %d and p = %d.\n', ...
    n, p);


% ============================================================
% 80-20 train/test split
% ============================================================

rng(seed);

perm = randperm(n);

n_train = floor(train_frac * n);

idx_train = perm(1:n_train);
idx_test  = perm(n_train+1:end);

X_train_raw = X_raw(idx_train, :);
X_test_raw  = X_raw(idx_test, :);

y_train_raw = y_raw(idx_train);
y_test_raw  = y_raw(idx_test);


% ============================================================
% Standardize X using TRAINING statistics only
% ============================================================

muX = mean(X_train_raw, 1);
sdX = std(X_train_raw, 0, 1);

% Protect zero-variance predictors
sdX(~isfinite(sdX) | sdX == 0) = 1;

X_train = (X_train_raw - muX) ./ sdX;
X_test  = (X_test_raw  - muX) ./ sdX;


% ============================================================
% Standardize y using TRAINING statistics only
% ============================================================

muY = mean(y_train_raw);
sdY = std(y_train_raw);

if sdY == 0
    error('Standard deviation of training y is zero.');
end

y_train = (y_train_raw - muY) / sdY;
y_test  = (y_test_raw  - muY) / sdY;


% ============================================================
% Convert to double
% ============================================================

X_train = double(X_train);
X_test  = double(X_test);

y_train = double(y_train);
y_test  = double(y_test);


% ============================================================
% SAVE SPLIT INFORMATION
% ============================================================

Tsplit = table();

Tsplit.obs_index = (1:n)';
Tsplit.in_train = ismember((1:n)', idx_train);
Tsplit.in_test  = ismember((1:n)', idx_test);

outfile_split = fullfile( ...
    outdir, ...
    'Output_8020_SIS_split.csv');

writetable(Tsplit, outfile_split);


% ============================================================
% STEP 1: SURE INDEPENDENCE SCREENING
%
% IMPORTANT:
% Use n_train, not total n.
% ============================================================

fprintf('\n===============================================\n');
fprintf('Running SIS screening...\n');
fprintf('===============================================\n');

tic;

d_sis = min(p, floor(n_train / log(n_train)));

% X_train and y_train have already been standardized.
%
% Therefore ranking by |X_j' y| is equivalent to ranking
% by absolute marginal correlation.
sis_score = abs(X_train' * y_train);

sis_score(~isfinite(sis_score)) = -Inf;

% More efficient than sorting all ~866K predictors
[~, sis_idx] = maxk(sis_score, d_sis);

X_train_sis = X_train(:, sis_idx);
X_test_sis  = X_test(:, sis_idx);

sis_time = toc;

fprintf('SIS retained %d of %d predictors.\n', d_sis, p);
fprintf('SIS time: %.3f seconds\n', sis_time);


% ============================================================
% Save common SIS screening information
% ============================================================

Tsis = table();

Tsis.feature_index = sis_idx(:);
Tsis.sis_score = sis_score(sis_idx);

% Sort saved file by SIS score
Tsis = sortrows(Tsis, 'sis_score', 'descend');

outfile_sis = fullfile( ...
    outdir, ...
    'Output_8020_SIS_screened_features.csv');

writetable(Tsis, outfile_sis);


% ############################################################
% ############################################################
%
%                   SIS-LASSO
%
% ############################################################
% ############################################################

fprintf('\n===============================================\n');
fprintf('Running SIS-LASSO...\n');
fprintf('===============================================\n');

rng(seed);

tic;

[B_lasso, FitInfo_lasso] = lasso( ...
    X_train_sis, ...
    y_train, ...
    'Alpha', 1, ...
    'CV', nfolds, ...
    'Standardize', true);

lasso_time = toc;

idx_lasso = FitInfo_lasso.IndexMinMSE;

beta_lasso_reduced = B_lasso(:, idx_lasso);
intercept_lasso = FitInfo_lasso.Intercept(idx_lasso);

lambda_lasso = FitInfo_lasso.Lambda(idx_lasso);

% Full p-dimensional coefficient vector
beta_hat_lasso = zeros(p,1);
beta_hat_lasso(sis_idx) = beta_lasso_reduced;

% ------------------------------------------------------------
% Test prediction
% ------------------------------------------------------------

yhat_lasso_std = ...
    intercept_lasso + X_test_sis * beta_lasso_reduced;

yhat_lasso_raw = ...
    muY + sdY * yhat_lasso_std;

resid_lasso_std = y_test - yhat_lasso_std;
resid_lasso_raw = y_test_raw - yhat_lasso_raw;

% ------------------------------------------------------------
% Metrics
% ------------------------------------------------------------

rmse_lasso = sqrt(mean( ...
    (y_test_raw - yhat_lasso_raw).^2));

mae_lasso = mean(abs( ...
    y_test_raw - yhat_lasso_raw));

if std(y_test_raw) > 0 && std(yhat_lasso_raw) > 0
    corr_lasso = corr(y_test_raw, yhat_lasso_raw);
else
    corr_lasso = NaN;
end

sst = sum((y_test_raw - mean(y_test_raw)).^2);
sse = sum((y_test_raw - yhat_lasso_raw).^2);

if sst > 0
    r2_lasso = 1 - sse/sst;
else
    r2_lasso = NaN;
end

n_selected_lasso = sum( ...
    abs(beta_hat_lasso) > sel_thresh);

total_time_lasso = sis_time + lasso_time;


% ------------------------------------------------------------
% Save LASSO predictions
% ------------------------------------------------------------

Tpred_lasso = table();

Tpred_lasso.obs_index = idx_test(:);
Tpred_lasso.y_true_raw = y_test_raw(:);
Tpred_lasso.y_pred_raw = yhat_lasso_raw(:);
Tpred_lasso.resid_raw = resid_lasso_raw(:);

Tpred_lasso.y_true_std = y_test(:);
Tpred_lasso.y_pred_std = yhat_lasso_std(:);
Tpred_lasso.resid_std = resid_lasso_std(:);

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_LASSO_test_predictions.csv');

writetable(Tpred_lasso, outfile);


% ------------------------------------------------------------
% Save LASSO metrics
% ------------------------------------------------------------

Tmetrics_lasso = table();

Tmetrics_lasso.n = n;
Tmetrics_lasso.p = p;

Tmetrics_lasso.n_train = n_train;
Tmetrics_lasso.n_test = n - n_train;

Tmetrics_lasso.d_sis = d_sis;

Tmetrics_lasso.rmse_test = rmse_lasso;
Tmetrics_lasso.mae_test = mae_lasso;
Tmetrics_lasso.corr_test = corr_lasso;
Tmetrics_lasso.r2_test = r2_lasso;

Tmetrics_lasso.n_selected = n_selected_lasso;

Tmetrics_lasso.lambda = lambda_lasso;

Tmetrics_lasso.sis_time_sec = sis_time;
Tmetrics_lasso.fit_time_sec = lasso_time;
Tmetrics_lasso.comp_time_sec = total_time_lasso;

Tmetrics_lasso.sel_thresh = sel_thresh;

Tmetrics_lasso.y_train_mean = muY;
Tmetrics_lasso.y_train_sd = sdY;

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_LASSO_metrics.csv');

writetable(Tmetrics_lasso, outfile);


% ------------------------------------------------------------
% Save LASSO screened coefficients
% ------------------------------------------------------------

Tcoef_lasso = table();

Tcoef_lasso.feature_index = sis_idx(:);
Tcoef_lasso.sis_score = sis_score(sis_idx);

Tcoef_lasso.beta_std = beta_lasso_reduced(:);

Tcoef_lasso.selected = ...
    abs(beta_lasso_reduced(:)) > sel_thresh;

Tcoef_lasso = sortrows( ...
    Tcoef_lasso, ...
    'sis_score', ...
    'descend');

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_LASSO_coefficients.csv');

writetable(Tcoef_lasso, outfile);


% ############################################################
% ############################################################
%
%                   SIS-SCAD
%
% ############################################################
% ############################################################

fprintf('\n===============================================\n');
fprintf('Running SIS-SCAD...\n');
fprintf('===============================================\n');

rng(seed);

tic;

[beta_scad_reduced, ...
 intercept_scad, ...
 lambda_scad] = fit_scad_cv( ...
    X_train_sis, ...
    y_train, ...
    nfolds, ...
    scad_gamma, ...
    seed, ...
    scad_num_lambda, ...
    scad_lambda_ratio, ...
    scad_max_iter, ...
    scad_tol);

scad_time = toc;

% Full p-dimensional coefficient vector
beta_hat_scad = zeros(p,1);
beta_hat_scad(sis_idx) = beta_scad_reduced;


% ------------------------------------------------------------
% Test prediction
% ------------------------------------------------------------

yhat_scad_std = ...
    intercept_scad + X_test_sis * beta_scad_reduced;

yhat_scad_raw = ...
    muY + sdY * yhat_scad_std;

resid_scad_std = y_test - yhat_scad_std;
resid_scad_raw = y_test_raw - yhat_scad_raw;


% ------------------------------------------------------------
% Metrics
% ------------------------------------------------------------

rmse_scad = sqrt(mean( ...
    (y_test_raw - yhat_scad_raw).^2));

mae_scad = mean(abs( ...
    y_test_raw - yhat_scad_raw));

if std(y_test_raw) > 0 && std(yhat_scad_raw) > 0
    corr_scad = corr(y_test_raw, yhat_scad_raw);
else
    corr_scad = NaN;
end

sst = sum((y_test_raw - mean(y_test_raw)).^2);
sse = sum((y_test_raw - yhat_scad_raw).^2);

if sst > 0
    r2_scad = 1 - sse/sst;
else
    r2_scad = NaN;
end

n_selected_scad = sum( ...
    abs(beta_hat_scad) > sel_thresh);

total_time_scad = sis_time + scad_time;


% ------------------------------------------------------------
% Save SCAD predictions
% ------------------------------------------------------------

Tpred_scad = table();

Tpred_scad.obs_index = idx_test(:);
Tpred_scad.y_true_raw = y_test_raw(:);
Tpred_scad.y_pred_raw = yhat_scad_raw(:);
Tpred_scad.resid_raw = resid_scad_raw(:);

Tpred_scad.y_true_std = y_test(:);
Tpred_scad.y_pred_std = yhat_scad_std(:);
Tpred_scad.resid_std = resid_scad_std(:);

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_SCAD_test_predictions.csv');

writetable(Tpred_scad, outfile);


% ------------------------------------------------------------
% Save SCAD metrics
% ------------------------------------------------------------

Tmetrics_scad = table();

Tmetrics_scad.n = n;
Tmetrics_scad.p = p;

Tmetrics_scad.n_train = n_train;
Tmetrics_scad.n_test = n - n_train;

Tmetrics_scad.d_sis = d_sis;

Tmetrics_scad.rmse_test = rmse_scad;
Tmetrics_scad.mae_test = mae_scad;
Tmetrics_scad.corr_test = corr_scad;
Tmetrics_scad.r2_test = r2_scad;

Tmetrics_scad.n_selected = n_selected_scad;

Tmetrics_scad.lambda = lambda_scad;
Tmetrics_scad.gamma_SCAD = scad_gamma;

Tmetrics_scad.sis_time_sec = sis_time;
Tmetrics_scad.fit_time_sec = scad_time;
Tmetrics_scad.comp_time_sec = total_time_scad;

Tmetrics_scad.sel_thresh = sel_thresh;

Tmetrics_scad.y_train_mean = muY;
Tmetrics_scad.y_train_sd = sdY;

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_SCAD_metrics.csv');

writetable(Tmetrics_scad, outfile);


% ------------------------------------------------------------
% Save SCAD screened coefficients
% ------------------------------------------------------------

Tcoef_scad = table();

Tcoef_scad.feature_index = sis_idx(:);
Tcoef_scad.sis_score = sis_score(sis_idx);

Tcoef_scad.beta_std = beta_scad_reduced(:);

Tcoef_scad.selected = ...
    abs(beta_scad_reduced(:)) > sel_thresh;

Tcoef_scad = sortrows( ...
    Tcoef_scad, ...
    'sis_score', ...
    'descend');

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_SCAD_coefficients.csv');

writetable(Tcoef_scad, outfile);


% ############################################################
% ############################################################
%
%                   SIS-RHS
%
% ############################################################
% ############################################################

fprintf('\n===============================================\n');
fprintf('Running SIS-RHS...\n');
fprintf('===============================================\n');

opts_rhs.seed = seed;

tic;

out_rhs = BUGSRHS_mcmc( ...
    X_train_sis, ...
    y_train, ...
    opts_rhs);

rhs_time = toc;


% ------------------------------------------------------------
% Posterior beta estimate
% ------------------------------------------------------------

switch lower(use_beta_summary)

    case 'mean'

        beta_rhs_reduced = out_rhs.beta_mean;

    case 'median'

        beta_rhs_reduced = out_rhs.beta_median;

    otherwise

        error(['use_beta_summary must be ' ...
            '''mean'' or ''median''.']);

end

beta_rhs_reduced = beta_rhs_reduced(:);


% Full p-dimensional coefficient vector
beta_hat_rhs = zeros(p,1);
beta_hat_rhs(sis_idx) = beta_rhs_reduced;


% ------------------------------------------------------------
% Test prediction
%
% RHS is fitted to centered / standardized data, so its
% standardized-scale intercept is zero.
% ------------------------------------------------------------

yhat_rhs_std = ...
    X_test_sis * beta_rhs_reduced;

yhat_rhs_raw = ...
    muY + sdY * yhat_rhs_std;

resid_rhs_std = y_test - yhat_rhs_std;
resid_rhs_raw = y_test_raw - yhat_rhs_raw;


% ------------------------------------------------------------
% Metrics
% ------------------------------------------------------------

rmse_rhs = sqrt(mean( ...
    (y_test_raw - yhat_rhs_raw).^2));

mae_rhs = mean(abs( ...
    y_test_raw - yhat_rhs_raw));

if std(y_test_raw) > 0 && std(yhat_rhs_raw) > 0
    corr_rhs = corr(y_test_raw, yhat_rhs_raw);
else
    corr_rhs = NaN;
end

sst = sum((y_test_raw - mean(y_test_raw)).^2);
sse = sum((y_test_raw - yhat_rhs_raw).^2);

if sst > 0
    r2_rhs = 1 - sse/sst;
else
    r2_rhs = NaN;
end

n_selected_rhs = sum( ...
    abs(beta_hat_rhs) > sel_thresh);

total_time_rhs = sis_time + rhs_time;


% ------------------------------------------------------------
% Save RHS predictions
% ------------------------------------------------------------

Tpred_rhs = table();

Tpred_rhs.obs_index = idx_test(:);
Tpred_rhs.y_true_raw = y_test_raw(:);
Tpred_rhs.y_pred_raw = yhat_rhs_raw(:);
Tpred_rhs.resid_raw = resid_rhs_raw(:);

Tpred_rhs.y_true_std = y_test(:);
Tpred_rhs.y_pred_std = yhat_rhs_std(:);
Tpred_rhs.resid_std = resid_rhs_std(:);

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_RHS_test_predictions.csv');

writetable(Tpred_rhs, outfile);


% ------------------------------------------------------------
% Save RHS metrics
% ------------------------------------------------------------

Tmetrics_rhs = table();

Tmetrics_rhs.n = n;
Tmetrics_rhs.p = p;

Tmetrics_rhs.n_train = n_train;
Tmetrics_rhs.n_test = n - n_train;

Tmetrics_rhs.d_sis = d_sis;

Tmetrics_rhs.beta_summary = ...
    string(use_beta_summary);

Tmetrics_rhs.rmse_test = rmse_rhs;
Tmetrics_rhs.mae_test = mae_rhs;
Tmetrics_rhs.corr_test = corr_rhs;
Tmetrics_rhs.r2_test = r2_rhs;

Tmetrics_rhs.n_selected = n_selected_rhs;

Tmetrics_rhs.sis_time_sec = sis_time;
Tmetrics_rhs.fit_time_sec = rhs_time;
Tmetrics_rhs.comp_time_sec = total_time_rhs;

Tmetrics_rhs.sel_thresh = sel_thresh;

Tmetrics_rhs.y_train_mean = muY;
Tmetrics_rhs.y_train_sd = sdY;

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_RHS_metrics.csv');

writetable(Tmetrics_rhs, outfile);


% ------------------------------------------------------------
% Save RHS screened coefficients
% ------------------------------------------------------------

Tcoef_rhs = table();

Tcoef_rhs.feature_index = sis_idx(:);
Tcoef_rhs.sis_score = sis_score(sis_idx);

Tcoef_rhs.beta_std = beta_rhs_reduced(:);

Tcoef_rhs.selected = ...
    abs(beta_rhs_reduced(:)) > sel_thresh;

Tcoef_rhs = sortrows( ...
    Tcoef_rhs, ...
    'sis_score', ...
    'descend');

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_RHS_coefficients.csv');

writetable(Tcoef_rhs, outfile);


% ============================================================
% RHS hyperparameter draws
% ============================================================

Thyper_rhs = table();

if isfield(out_rhs, 'sigma2_draws')
    Thyper_rhs.sigma2 = out_rhs.sigma2_draws(:);
end

if isfield(out_rhs, 'tau_draws')
    Thyper_rhs.tau = out_rhs.tau_draws(:);
end

if isfield(out_rhs, 'c2_draws')
    Thyper_rhs.c2 = out_rhs.c2_draws(:);
end

if isfield(out_rhs, 'eta_draws')
    Thyper_rhs.eta = out_rhs.eta_draws(:);
end

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_RHS_hyper_draws.csv');

writetable(Thyper_rhs, outfile);


% ============================================================
% RHS hyperparameter summaries
% ============================================================

param_names = {};
post_mean = [];
post_median = [];
post_sd = [];
ci_2p5 = [];
ci_97p5 = [];


if isfield(out_rhs, 'sigma2_draws')

    x = out_rhs.sigma2_draws(:);

    param_names{end+1,1} = 'sigma2';

    post_mean(end+1,1) = mean(x);
    post_median(end+1,1) = median(x);
    post_sd(end+1,1) = std(x);

    ci_2p5(end+1,1) = prctile(x,2.5);
    ci_97p5(end+1,1) = prctile(x,97.5);

end


if isfield(out_rhs, 'tau_draws')

    x = out_rhs.tau_draws(:);

    param_names{end+1,1} = 'tau';

    post_mean(end+1,1) = mean(x);
    post_median(end+1,1) = median(x);
    post_sd(end+1,1) = std(x);

    ci_2p5(end+1,1) = prctile(x,2.5);
    ci_97p5(end+1,1) = prctile(x,97.5);

end


if isfield(out_rhs, 'c2_draws')

    x = out_rhs.c2_draws(:);

    param_names{end+1,1} = 'c2';

    post_mean(end+1,1) = mean(x);
    post_median(end+1,1) = median(x);
    post_sd(end+1,1) = std(x);

    ci_2p5(end+1,1) = prctile(x,2.5);
    ci_97p5(end+1,1) = prctile(x,97.5);

end


if isfield(out_rhs, 'eta_draws')

    x = out_rhs.eta_draws(:);

    param_names{end+1,1} = 'eta';

    post_mean(end+1,1) = mean(x);
    post_median(end+1,1) = median(x);
    post_sd(end+1,1) = std(x);

    ci_2p5(end+1,1) = prctile(x,2.5);
    ci_97p5(end+1,1) = prctile(x,97.5);

end


Thyper_summary_rhs = table( ...
    param_names, ...
    post_mean, ...
    post_median, ...
    post_sd, ...
    ci_2p5, ...
    ci_97p5, ...
    'VariableNames', ...
    {'parameter','mean','median','sd','ci_2p5','ci_97p5'});

outfile = fullfile( ...
    outdir, ...
    'Output_8020_SIS_RHS_hyper_summary.csv');

writetable(Thyper_summary_rhs, outfile);


% ============================================================
% Console summary
% ============================================================

disp(' ');
disp('====================================================');
disp('80-20 SIS ANALYSIS COMPLETED');
disp('====================================================');

fprintf('n       = %d\n', n);
fprintf('p       = %d\n', p);
fprintf('n_train = %d\n', n_train);
fprintf('n_test  = %d\n', n-n_train);
fprintf('d_SIS   = %d\n', d_sis);

disp(' ');
disp('SIS-LASSO:');
disp(Tmetrics_lasso);

disp(' ');
disp('SIS-SCAD:');
disp(Tmetrics_scad);

disp(' ');
disp('SIS-RHS:');
disp(Tmetrics_rhs);

disp(' ');
fprintf('Saved files to folder: %s\n', outdir);
disp('====================================================');


% ============================================================
% ============================================================
% LOCAL FUNCTIONS FOR SCAD
% ============================================================
% ============================================================


function [beta_best, intercept_best, lambda_best] = ...
    fit_scad_cv(X, y, K, gamma_scad, seed, ...
                num_lambda, lambda_ratio, max_iter, tol)

% Gaussian SCAD with K-fold cross-validation.

rng(seed);

[n, ~] = size(X);
y = y(:);


% ------------------------------------------------------------
% Common lambda grid based on complete training data
% ------------------------------------------------------------

muX = mean(X,1);

Xc = X - muX;

xscale = sqrt(sum(Xc.^2,1) / n);

xscale(~isfinite(xscale) | xscale <= 0) = 1;

Z = Xc ./ xscale;

yc = y - mean(y);

lambda_max = max(abs(Z' * yc)) / n;

if ~isfinite(lambda_max) || lambda_max <= 0

    beta_best = zeros(size(X,2),1);
    intercept_best = mean(y);
    lambda_best = 0;

    return;

end

lambda_min = lambda_max * lambda_ratio;

lambda_seq = exp(linspace( ...
    log(lambda_max), ...
    log(lambda_min), ...
    num_lambda));


% ------------------------------------------------------------
% Balanced random folds
% ------------------------------------------------------------

perm = randperm(n);

fold_id = zeros(n,1);

for ii = 1:n
    fold_id(perm(ii)) = mod(ii-1,K) + 1;
end


% ------------------------------------------------------------
% Cross-validation
% ------------------------------------------------------------

cv_mse = zeros(K, num_lambda);

for k = 1:K

    train_idx = fold_id ~= k;
    test_idx  = fold_id == k;

    Xtr = X(train_idx,:);
    ytr = y(train_idx);

    Xva = X(test_idx,:);
    yva = y(test_idx);

    [B_path, intercept_path] = scad_path_fit( ...
        Xtr, ...
        ytr, ...
        lambda_seq, ...
        gamma_scad, ...
        max_iter, ...
        tol);

    pred = Xva * B_path;

    pred = pred + repmat( ...
        intercept_path, ...
        size(Xva,1), ...
        1);

    err = yva - pred;

    cv_mse(k,:) = mean(err.^2,1);

end


% ------------------------------------------------------------
% Minimum CV error
% ------------------------------------------------------------

mean_cv_mse = mean(cv_mse,1);

[~, best_idx] = min(mean_cv_mse);

lambda_best = lambda_seq(best_idx);


% ------------------------------------------------------------
% Refit complete training data
% ------------------------------------------------------------

[B_full, intercept_full] = scad_path_fit( ...
    X, ...
    y, ...
    lambda_seq, ...
    gamma_scad, ...
    max_iter, ...
    tol);

beta_best = B_full(:,best_idx);
intercept_best = intercept_full(best_idx);

end


% ============================================================

function [B_original, intercepts] = ...
    scad_path_fit(X, y, lambda_seq, gamma_scad, max_iter, tol)

[n, pp] = size(X);

y = y(:);


% ------------------------------------------------------------
% Standardize predictors
% ------------------------------------------------------------

muX = mean(X,1);

Xc = X - muX;

xscale = sqrt(sum(Xc.^2,1) / n);

invalid_scale = ...
    ~isfinite(xscale) | xscale <= 0;

xscale(invalid_scale) = 1;

Z = Xc ./ xscale;


% ------------------------------------------------------------
% Center response
% ------------------------------------------------------------

muy = mean(y);
yc = y - muy;


% ------------------------------------------------------------
% Initialization
% ------------------------------------------------------------

num_lambda = length(lambda_seq);

B_standardized = zeros(pp,num_lambda);

beta = zeros(pp,1);

resid = yc;


% ------------------------------------------------------------
% Lambda path
% ------------------------------------------------------------

for ll = 1:num_lambda

    lambda = lambda_seq(ll);

    for iter = 1:max_iter

        max_change = 0;

        for j = 1:pp

            beta_old = beta(j);

            % Add previous contribution back
            resid = resid + Z(:,j) * beta_old;

            % Standardized coordinate score
            z_j = (Z(:,j)' * resid) / n;

            beta_new = scad_threshold( ...
                z_j, ...
                lambda, ...
                gamma_scad);

            beta(j) = beta_new;

            % Remove new contribution
            resid = resid - Z(:,j) * beta_new;

            change = abs(beta_new - beta_old);

            if change > max_change
                max_change = change;
            end

        end

        if max_change < tol
            break;
        end

    end

    B_standardized(:,ll) = beta;

end


% ------------------------------------------------------------
% Return to original predictor scale
% ------------------------------------------------------------

B_original = ...
    B_standardized ./ xscale(:);

B_original(invalid_scale,:) = 0;

intercepts = ...
    muy - muX * B_original;

end


% ============================================================

function beta = scad_threshold(z, lambda, a)

az = abs(z);

if az <= lambda

    beta = 0;

elseif az <= 2*lambda

    beta = sign(z) * (az-lambda);

elseif az <= a*lambda

    beta = ...
        ((a-1)*z - sign(z)*a*lambda) / ...
        (a-2);

else

    beta = z;

end

end