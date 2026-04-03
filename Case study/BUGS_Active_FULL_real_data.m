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
%p = keep_p_upto;

%% =========================
%  BUGS options
%  =========================
opts.niter = 5000;
opts.burnin = 1000;
opts.thin = 5;
opts.display_every = 2;
opts.verbose = true;

opts.tau0 = 0.0005;           % smaller => stronger global shrinkage
opts.sigma_eta = 1e-6;       % guidance useful => put higher
opts.Cz = 2;
opts.a_c = 10;               % E(c^2)=b_c/(a_c-1), small c^2 => more shrinkage
opts.b_c = 100;
opts.a_sigma = 2;
opts.b_sigma = 2;

%% =========================
%  BUGS-Active specific options
%  =========================
opts.save_full_beta = true;          % needed here for full posterior summaries / CI
opts.max_p_full_save = p;            % must be >= p
opts.save_top_beta = min(p, 2000);    % irrelevant, post summary uses full draw

% active-set controls
opts.active_beta_thresh = 0.0005;
opts.active_cap = min(p, 20000);
opts.min_active = min(p, 10000);
opts.guidance_keep = min(p, 10000);
opts.lambda2_inactive = 1e-3;

%% =========================
%  User settings
%  =========================
seed = 1;
rng(seed);

use_beta_summary = 'mean';   % 'mean' or 'median'

% Phase 2 output controls
sel_thresh = 0.01;
topK = 200;   % for top-CpG tables / forest plot prep

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
%  Standardize X using full data
%  =========================
muX = mean(X_raw, 1);
sdX = std(X_raw, 0, 1);

% protect against zero-variance columns
sdX(sdX == 0) = 1;

X = (X_raw - muX) ./ sdX;

%% =========================
%  Standardize y using full data
%  =========================
muY = mean(y_raw);
sdY = std(y_raw);

if sdY == 0
    error('Standard deviation of y is zero.');
end

y = (y_raw - muY) / sdY;

%% =========================
%  Converting to double for optimized BLAS
%  =========================
X = double(X);
y = double(y);

%% =========================
%  Fit BUGS_Active_mcmc on full data
%  =========================
tic;
out = BUGS_Active_mcmc(X, y, opts);
comp_time = toc;

%% =========================
%  Choose beta estimate
%  =========================
switch lower(use_beta_summary)
    case 'mean'
        beta_hat = out.beta_mean;
        ranking_metric = abs(out.beta_mean(:));
    case 'median'
        beta_hat = out.beta_median;
        ranking_metric = abs(out.beta_median(:));
    otherwise
        error('use_beta_summary must be ''mean'' or ''median''.');
end

beta_hat = beta_hat(:);

if length(beta_hat) ~= p
    error('beta_hat has length %d, expected %d.', length(beta_hat), p);
end

%% =========================
%  Beta posterior summaries
%  =========================
if ~isfield(out, 'beta_draws') || isempty(out.beta_draws)
    error('Full beta draws not found in output. Phase 2 requires opts.save_full_beta = true.');
end

beta_draws = out.beta_draws;                 % nkeep x p
beta_mean = out.beta_mean(:);
beta_median = out.beta_median(:);
beta_sd = std(beta_draws, 0, 1)';
beta_ci_2p5 = prctile(beta_draws, 2.5, 1)';
beta_ci_97p5 = prctile(beta_draws, 97.5, 1)';
post_prob_abs_gt_thresh = mean(abs(beta_draws) > sel_thresh, 1)';

selected_mean_thresh = abs(beta_mean) > sel_thresh;
selected_median_thresh = abs(beta_median) > sel_thresh;
selected_ci_excludes_zero = (beta_ci_2p5 > 0) | (beta_ci_97p5 < 0);

n_selected_mean = sum(selected_mean_thresh);
n_selected_median = sum(selected_median_thresh);
n_selected_ci = sum(selected_ci_excludes_zero);

%% =========================
%  Raw-scale beta summaries
%  y_std = X_std * beta  => raw-scale coefficient on X_raw:
%  beta_raw = sdY * beta_std ./ sdX
%  intercept_raw = muY - sum(muX .* beta_raw)
%  We save raw-scale slopes too.
%  =========================
scale_factor = (sdY ./ sdX(:));
beta_mean_raw = beta_mean .* scale_factor;
beta_median_raw = beta_median .* scale_factor;
beta_ci_2p5_raw = beta_ci_2p5 .* scale_factor;
beta_ci_97p5_raw = beta_ci_97p5 .* scale_factor;

intercept_mean_raw = muY - sum(muX(:) .* beta_mean_raw);
intercept_median_raw = muY - sum(muX(:) .* beta_median_raw);

%% =========================
%  Top-K CpGs for paper plots/tables
%  =========================
topK = min(topK, p);
[~, ord_top] = sort(ranking_metric, 'descend');
idx_top = ord_top(1:topK);

Ttop = table();
Ttop.rank = (1:topK)';
Ttop.index = idx_top(:);
Ttop.beta_mean = beta_mean(idx_top);
Ttop.beta_median = beta_median(idx_top);
Ttop.beta_sd = beta_sd(idx_top);
Ttop.ci_2p5 = beta_ci_2p5(idx_top);
Ttop.ci_97p5 = beta_ci_97p5(idx_top);
Ttop.beta_mean_raw = beta_mean_raw(idx_top);
Ttop.beta_median_raw = beta_median_raw(idx_top);
Ttop.ci_2p5_raw = beta_ci_2p5_raw(idx_top);
Ttop.ci_97p5_raw = beta_ci_97p5_raw(idx_top);
Ttop.abs_beta_mean = abs(beta_mean(idx_top));
Ttop.abs_beta_median = abs(beta_median(idx_top));
Ttop.post_prob_abs_gt_thresh = post_prob_abs_gt_thresh(idx_top);
Ttop.selected_mean_thresh = selected_mean_thresh(idx_top);
Ttop.selected_median_thresh = selected_median_thresh(idx_top);
Ttop.selected_ci_excludes_zero = selected_ci_excludes_zero(idx_top);
Ttop.sign_beta_mean = sign(beta_mean(idx_top));
Ttop.sign_beta_median = sign(beta_median(idx_top));

outfile_top = fullfile(outdir, 'Output_FULL_topK_beta_summary.csv');
writetable(Ttop, outfile_top);

%% =========================
%  Full beta summary table for all CpGs
%  =========================
Tbeta = table();
Tbeta.index = (1:p)';
Tbeta.beta_mean = beta_mean;
Tbeta.beta_median = beta_median;
Tbeta.beta_sd = beta_sd;
Tbeta.ci_2p5 = beta_ci_2p5;
Tbeta.ci_97p5 = beta_ci_97p5;
Tbeta.beta_mean_raw = beta_mean_raw;
Tbeta.beta_median_raw = beta_median_raw;
Tbeta.ci_2p5_raw = beta_ci_2p5_raw;
Tbeta.ci_97p5_raw = beta_ci_97p5_raw;
Tbeta.abs_beta_mean = abs(beta_mean);
Tbeta.abs_beta_median = abs(beta_median);
Tbeta.post_prob_abs_gt_thresh = post_prob_abs_gt_thresh;
Tbeta.selected_mean_thresh = selected_mean_thresh;
Tbeta.selected_median_thresh = selected_median_thresh;
Tbeta.selected_ci_excludes_zero = selected_ci_excludes_zero;

outfile_beta = fullfile(outdir, 'Output_FULL_beta_summary.csv');
writetable(Tbeta, outfile_beta);

%% =========================
%  Save top-K posterior draws only (useful for forest plot / diagnostics)
%  =========================
Ttop_draws = array2table(beta_draws(:, idx_top));
top_draw_names = strcat("beta_draw_top_rank_", string(1:topK));
Ttop_draws.Properties.VariableNames = cellstr(top_draw_names);

outfile_top_draws = fullfile(outdir, 'Output_FULL_topK_beta_draws.csv');
writetable(Ttop_draws, outfile_top_draws);

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

outfile_hdraws = fullfile(outdir, 'Output_FULL_hyper_draws.csv');
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

outfile_hsummary = fullfile(outdir, 'Output_FULL_hyper_summary.csv');
writetable(Thyper_summary, outfile_hsummary);

%% =========================
%  Save metadata / run summary
%  =========================
Tmeta = table();
Tmeta.n = n;
Tmeta.p = p;
Tmeta.seed = seed;
Tmeta.keep_n_upto = string(num2str(n_keep));
Tmeta.keep_p_upto = string(num2str(p_keep));
Tmeta.beta_summary = string(use_beta_summary);
Tmeta.sel_thresh = sel_thresh;
Tmeta.topK = topK;
Tmeta.comp_time_sec = comp_time;
Tmeta.y_mean_raw = muY;
Tmeta.y_sd_raw = sdY;
Tmeta.intercept_mean_raw = intercept_mean_raw;
Tmeta.intercept_median_raw = intercept_median_raw;
Tmeta.n_selected_mean = n_selected_mean;
Tmeta.n_selected_median = n_selected_median;
Tmeta.n_selected_ci = n_selected_ci;
Tmeta.niter = opts.niter;
Tmeta.burnin = opts.burnin;
Tmeta.thin = opts.thin;
Tmeta.tau0 = opts.tau0;
Tmeta.sigma_eta = opts.sigma_eta;
Tmeta.Cz = opts.Cz;
Tmeta.a_c = opts.a_c;
Tmeta.b_c = opts.b_c;
Tmeta.a_sigma = opts.a_sigma;
Tmeta.b_sigma = opts.b_sigma;
Tmeta.save_full_beta = opts.save_full_beta;
Tmeta.max_p_full_save = opts.max_p_full_save;
Tmeta.save_top_beta = opts.save_top_beta;
Tmeta.active_beta_thresh = opts.active_beta_thresh;
Tmeta.active_cap = opts.active_cap;
Tmeta.min_active = opts.min_active;
Tmeta.guidance_keep = opts.guidance_keep;
Tmeta.lambda2_inactive = opts.lambda2_inactive;

outfile_meta = fullfile(outdir, 'Output_FULL_metadata.csv');
writetable(Tmeta, outfile_meta);

%% =========================
%  Optional: save a MAT file with the main objects
%  Uncomment if you want an archive for later without rerunning
%  =========================
% save(fullfile(outdir, 'Output_FULL_main_objects.mat'), ...
%     'out', 'idx_top', 'muY', 'sdY', 'muX', 'sdX', ...
%     'beta_mean', 'beta_median', 'beta_ci_2p5', 'beta_ci_97p5', ...
%     '-v7.3');

%% =========================
%  Console output
%  =========================
disp('====================================================');
disp('BUGS_Active_mcmc FULL-DATA run completed.');

disp(Tmeta);
disp(Thyper_summary);

fprintf('Saved files to folder: %s\n', outdir);
disp('Main Phase 2 outputs saved as:');
disp('  Output_FULL_beta_summary.csv');
disp('  Output_FULL_topK_beta_summary.csv');
disp('  Output_FULL_topK_beta_draws.csv');
disp('  Output_FULL_hyper_draws.csv');
disp('  Output_FULL_hyper_summary.csv');
disp('  Output_FULL_metadata.csv');
disp('====================================================');