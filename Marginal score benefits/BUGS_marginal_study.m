clearvars;
addpath('./BUGS/');
addpath('./Data/');
addpath('./supp funs/');

% ============================================================
% Scenario settings
% ============================================================
n = 100;
p = 10;
rep = 1;

% ============================================================
% BUGS options
% ============================================================
opts.niter = 5000;
opts.burnin = 1000;
opts.thin = 5;
opts.display_every = 100;
opts.verbose = true;

opts.tau0 = 0.001;           % smaller => stronger global shrinkage
opts.sigma_eta = 1;        % guidance strength prior scale
opts.Cz = 4;
opts.a_c = 2;               % E(c^2)=b_c/(a_c-1)
opts.b_c = 8;
opts.a_sigma = 2;
opts.b_sigma = 2;

opts.eps_score = 1e-6;      
opts.slice_w = 0.5;
opts.slice_m = 100;
opts.seed = 1;

% IMPORTANT: save all lambda draws since p is small here
opts.save_lambda_max = p;

% Posterior summary choice for overall metric table
use_beta_summary = 'mean';   % 'mean' or 'median'
sel_thresh = 0.01;           % threshold for posterior exceedance / selection summaries

fprintf('performing data rep: %d | n: %d | p: %d\n', rep, n, p);

% Create Results folder if needed
if ~exist('Results', 'dir')
    mkdir('Results');
end

% ============================================================
% File names
% ============================================================
fileX = sprintf('Data_marginal_X_n_%d_p_%d_rep_%d.csv', n, p, rep);
fileY = sprintf('Data_marginal_Y_n_%d_p_%d_rep_%d.csv', n, p, rep);
fileB = sprintf('Data_marginal_beta_n_%d_p_%d.csv', n, p);

% ============================================================
% Read data
% ============================================================
X = readmatrix(fileX);
y_raw = readmatrix(fileY);
beta_true = readmatrix(fileB);

y_raw = y_raw(:);
beta_true = beta_true(:);

if size(X,1) ~= n || size(X,2) ~= p
    error('X has dimension %d x %d, expected %d x %d.', size(X,1), size(X,2), n, p);
end
if length(y_raw) ~= n
    error('y has length %d, expected %d.', length(y_raw), n);
end
if length(beta_true) ~= p
    error('beta_true has length %d, expected %d.', length(beta_true), p);
end

% ============================================================
% Standardize y and corresponding beta_true
% ============================================================
y_raw = y_raw - mean(y_raw);
sy = std(y_raw);

if sy == 0
    error('Standard deviation of y_raw is zero.');
end

y = y_raw / sy;
beta_true_std = beta_true / sy;

% ============================================================
% Compute marginal guidance statistics outside BUGS_mcmc
% (for saving and plotting)
% Keep this consistent with BUGS_mcmc
% ============================================================
s = abs(X' * y) / n;
ztilde = log(s + opts.eps_score);
z = (ztilde - mean(ztilde)) / std(ztilde);
z_clip = max(-opts.Cz, min(z, opts.Cz));

var_idx = (1:p)';
is_signal = abs(beta_true_std) > 0;

% ============================================================
% Save marginal score table
% Panel 1 of intro figure can use this directly
% ============================================================
T_scores = table( ...
    var_idx, beta_true, beta_true_std, is_signal, ...
    s, ztilde, z, z_clip, ...
    'VariableNames', { ...
        'index', 'beta_true_raw', 'beta_true_std', 'is_signal', ...
        'marginal_score', 'log_marginal_score', 'z_score', 'z_score_clipped' ...
    });

outfile_scores = sprintf('Results/Output_marginal_scores_n_%d_p_%d_rep_%d.csv', n, p, rep);
writetable(T_scores, outfile_scores);
fprintf('Saved marginal score table to: %s\n', outfile_scores);

% ============================================================
% Fit BUGS and record time
% ============================================================
tic;
out = BUGS_mcmc(X, y, opts);
comp_time = toc;

% Use BUGS-generated zstar for posterior calculations
if max(abs(out.zstar(:) - z_clip(:))) > 1e-10
    warning('External z_clip and out.zstar differ slightly. Using out.zstar from BUGS output.');
end
zstar = out.zstar(:);

% ============================================================
% Choose beta estimate for one-row overall performance metrics
% ============================================================
switch lower(use_beta_summary)
    case 'mean'
        beta_hat = out.beta_mean(:);
    case 'median'
        beta_hat = out.beta_median(:);
    otherwise
        error('use_beta_summary must be ''mean'' or ''median''.');
end

% ============================================================
% Compute overall metrics
% ============================================================
metrics = compute_metrics_from_beta_BUGS(X, beta_true_std, beta_hat, comp_time, sel_thresh);
T_metrics = struct2table(metrics);

outfile_metrics = sprintf('Results/Output_marginal_n_%d_p_%d_BUGS_rep_%d.csv', n, p, rep);
writetable(T_metrics, outfile_metrics);
disp(T_metrics);
fprintf('Saved overall metrics to: %s\n', outfile_metrics);

% ============================================================
% Extract posterior draws from current BUGS_mcmc output
% ============================================================
beta_draws = out.beta_draws;                 % nSave x p
tau_draws = out.tau_draws(:);               % stored as tau
c2_draws = out.c2_draws(:);                 % nSave x 1
eta_draws = out.eta_draws(:);               % nSave x 1
lambda_draws = out.lambda_draws_partial;    % nSave x p because save_lambda_max = p

[nSave, p_beta] = size(beta_draws);
[nSave_lambda, p_lambda] = size(lambda_draws);

if p_beta ~= p
    error('beta_draws has %d columns, expected p=%d.', p_beta, p);
end
if nSave_lambda ~= nSave || p_lambda ~= p
    error('lambda_draws_partial has size %d x %d, expected %d x %d.', ...
        nSave_lambda, p_lambda, nSave, p);
end
if length(tau_draws) ~= nSave || length(c2_draws) ~= nSave || length(eta_draws) ~= nSave
    error('tau_draws, c2_draws, eta_draws must all have length nSave.');
end

% ============================================================
% Reconstruct guided BUGS shrinkage quantities
% ============================================================
lambda2_draws = lambda_draws .^ 2;          % nSave x p
tau2_draws = tau_draws .^ 2;                % nSave x 1

zstar_row = reshape(zstar, 1, p);
guided_mult_draws = exp(eta_draws * zstar_row);                % nSave x p
base_scale_draws = lambda2_draws .* tau2_draws;                % nSave x p
guided_scale_draws = base_scale_draws .* guided_mult_draws;    % nSave x p

c2_mat = repmat(c2_draws, 1, p);

% Guided BUGS effective scale
kappa2_draws = (c2_mat .* guided_scale_draws) ./ (c2_mat + guided_scale_draws);
kappa2_draws = max(kappa2_draws, 1e-12);

% ============================================================
% Variable-wise posterior summaries
% Useful for intro panel 3 and section 2 summaries
% ============================================================
beta_mean = mean(beta_draws, 1)';
beta_median = median(beta_draws, 1)';
beta_sd = std(beta_draws, 0, 1)';

ci_2p5 = quantile(beta_draws, 0.025, 1)';
ci_97p5 = quantile(beta_draws, 0.975, 1)';

abs_beta_mean = mean(abs(beta_draws), 1)';
abs_beta_median = median(abs(beta_draws), 1)';

post_prob_abs_gt_thresh = mean(abs(beta_draws) > sel_thresh, 1)';

selected_mean_thresh = abs(beta_mean) > sel_thresh;
selected_median_thresh = abs(beta_median) > sel_thresh;
selected_ci_excludes_zero = (ci_2p5 > 0) | (ci_97p5 < 0);

sign_beta_mean = sign(beta_mean);
sign_beta_median = sign(beta_median);

lambda2_mean = mean(lambda2_draws, 1)';
lambda2_median = median(lambda2_draws, 1)';

guided_scale_mean = mean(guided_scale_draws, 1)';
guided_scale_median = median(guided_scale_draws, 1)';

kappa2_mean = mean(kappa2_draws, 1)';
kappa2_median = median(kappa2_draws, 1)';

T_varsum = table( ...
    var_idx, beta_true, beta_true_std, is_signal, ...
    s, ztilde, z, zstar, ...
    beta_mean, beta_median, beta_sd, ci_2p5, ci_97p5, ...
    abs_beta_mean, abs_beta_median, ...
    post_prob_abs_gt_thresh, ...
    selected_mean_thresh, selected_median_thresh, selected_ci_excludes_zero, ...
    sign_beta_mean, sign_beta_median, ...
    lambda2_mean, lambda2_median, ...
    guided_scale_mean, guided_scale_median, ...
    kappa2_mean, kappa2_median, ...
    'VariableNames', { ...
        'index', 'beta_true_raw', 'beta_true_std', 'is_signal', ...
        'marginal_score', 'log_marginal_score', 'z_score', 'z_score_clipped', ...
        'beta_mean', 'beta_median', 'beta_sd', 'ci_2p5', 'ci_97p5', ...
        'abs_beta_mean', 'abs_beta_median', ...
        'post_prob_abs_gt_thresh', ...
        'selected_mean_thresh', 'selected_median_thresh', 'selected_ci_excludes_zero', ...
        'sign_beta_mean', 'sign_beta_median', ...
        'lambda2_mean', 'lambda2_median', ...
        'guided_scale_mean', 'guided_scale_median', ...
        'kappa2_mean', 'kappa2_median' ...
    });

outfile_varsum = sprintf('Results/Output_BUGS_variable_summary_n_%d_p_%d_rep_%d.csv', n, p, rep);
writetable(T_varsum, outfile_varsum);
fprintf('Saved variable summary table to: %s\n', outfile_varsum);

% ============================================================
% Save long-format draws for boxplots / violin plots
% This is the main file for section-2 and intro figure plotting
% ============================================================
draw_id_col = repelem((1:nSave)', p, 1);
index_col = repmat((1:p)', nSave, 1);
is_signal_col = repmat(is_signal(:), nSave, 1);
marginal_score_col = repmat(s(:), nSave, 1);
zstar_col = repmat(zstar(:), nSave, 1);

beta_draw_col = beta_draws(:);
abs_beta_draw_col = abs(beta_draw_col);

lambda_draw_col = lambda_draws(:);
lambda2_draw_col = lambda2_draws(:);

tau_draw_col = repelem(tau_draws, p, 1);
tau2_draw_col = repelem(tau2_draws, p, 1);

c2_draw_col = repelem(c2_draws, p, 1);
eta_draw_col = repelem(eta_draws, p, 1);

base_scale_draw_col = base_scale_draws(:);
guided_scale_draw_col = guided_scale_draws(:);
kappa2_draw_col = kappa2_draws(:);

T_draws = table( ...
    draw_id_col, index_col, is_signal_col, marginal_score_col, zstar_col, ...
    beta_draw_col, abs_beta_draw_col, ...
    lambda_draw_col, lambda2_draw_col, ...
    tau_draw_col, tau2_draw_col, ...
    c2_draw_col, eta_draw_col, ...
    base_scale_draw_col, guided_scale_draw_col, ...
    kappa2_draw_col, ...
    'VariableNames', { ...
        'draw_id', 'index', 'is_signal', 'marginal_score', 'z_score_clipped', ...
        'beta_draw', 'abs_beta_draw', ...
        'lambda_draw', 'lambda2_draw', ...
        'tau_draw', 'tau2_draw', ...
        'c2_draw', 'eta_draw', ...
        'base_scale_draw', 'guided_scale_draw', ...
        'kappa2_draw' ...
    });

outfile_draws = sprintf('Results/Output_BUGS_shrinkage_draws_n_%d_p_%d_rep_%d.csv', n, p, rep);
writetable(T_draws, outfile_draws);
fprintf('Saved long-format draw table to: %s\n', outfile_draws);

% ============================================================
% Save hyperparameter posterior summary
% ============================================================
T_hyper = table( ...
    mean(tau_draws), median(tau_draws), quantile(tau_draws, 0.025), quantile(tau_draws, 0.975), ...
    mean(c2_draws), median(c2_draws), quantile(c2_draws, 0.025), quantile(c2_draws, 0.975), ...
    mean(eta_draws), median(eta_draws), quantile(eta_draws, 0.025), quantile(eta_draws, 0.975), ...
    'VariableNames', { ...
        'tau_mean', 'tau_median', 'tau_ci_2p5', 'tau_ci_97p5', ...
        'c2_mean', 'c2_median', 'c2_ci_2p5', 'c2_ci_97p5', ...
        'eta_mean', 'eta_median', 'eta_ci_2p5', 'eta_ci_97p5' ...
    });

outfile_hyper = sprintf('Results/Output_BUGS_hyper_summary_n_%d_p_%d_rep_%d.csv', n, p, rep);
writetable(T_hyper, outfile_hyper);
fprintf('Saved hyperparameter summary to: %s\n', outfile_hyper);

% ============================================================
% Save a compact section-2-specific summary file
% ============================================================
T_section2 = table( ...
    var_idx, is_signal, s, zstar, ...
    lambda2_mean, guided_scale_mean, kappa2_mean, ...
    post_prob_abs_gt_thresh, ...
    'VariableNames', { ...
        'index', 'is_signal', 'marginal_score', 'z_score_clipped', ...
        'lambda2_mean', 'guided_scale_mean', 'kappa2_mean', ...
        'post_prob_abs_gt_thresh' ...
    });

outfile_sec2 = sprintf('Results/Output_BUGS_section2_summary_n_%d_p_%d_rep_%d.csv', n, p, rep);
writetable(T_section2, outfile_sec2);
fprintf('Saved Section 2 summary table to: %s\n', outfile_sec2);

disp('All guided BUGS plot-ready outputs saved successfully.');