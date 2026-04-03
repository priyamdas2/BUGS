clearvars;
addpath('./BUGS/');
addpath('./Data/');
addpath('./supp funs/');

% ============================================================
% Scenario settings
% ============================================================
n = 200;
p = 1000;
rep = 1;
Active_cutoff = 0.0001;
guidance_keep_val = 20;

%% BUGS options

opts.niter = 2000;
opts.burnin = 500;
opts.thin = 5;
opts.display_each = 100;
opts.verbose = true;

opts.tau0 = 0.001;       % smaller => stronger global shrinkage
opts.sigma_eta = 0.002;  % guidance useful => put higher
opts.Cz = 2;
opts.a_c = 10;   % E(c^2) =b_c/(a_c-1), small c^2 => more shrinkage, small a_c  more heavy tail
opts.b_c = 100;
opts.a_sigma = 2;
opts.b_sigma = 2;

%% BUGS-Active specific options
opts.save_full_beta = true;      % needed for out2.beta_mean / out2.beta_median
opts.max_p_full_save = p;        % must be >= p

% active-set controls
opts.active_beta_thresh = Active_cutoff;
opts.active_cap = p;
% opts.min_active = min(p,10);
opts.guidance_keep = min(p,guidance_keep_val);
opts.lambda2_inactive = 1e-3;

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
fileX = sprintf('Data_cutoffEffect_X_n_%d_p_%d_rep_%d.csv', n, p, rep);
fileY = sprintf('Data_cutoffEffect_Y_n_%d_p_%d_rep_%d.csv', n, p, rep);
fileB = sprintf('Data_cutoffEffect_beta_n_%d_p_%d.csv', n, p);

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
% Fit BUGS-Active and record time
% ============================================================
tic;
out = BUGS_Active_mcmc(X, y, opts);
comp_time = toc;

% ============================================================
% User settings
% ============================================================
true_p = 5;                  % save posterior draws for first 5 true signals

% ============================================================
% Choose posterior point estimate for metrics
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
% Compute compact metrics
% ============================================================
metrics = compute_metrics_from_beta_BUGS(X, beta_true_std, beta_hat, comp_time, sel_thresh);

T_metrics = table( ...
    metrics.RMSE_beta, metrics.PredMSE, metrics.TPR, metrics.FPR, metrics.FDR, metrics.CompTime, ...
    'VariableNames', {'RMSE(beta)', 'MSE(y)', 'TPR', 'FPR', 'FDR', 'CompTime'} ...
);

outfile_metrics = sprintf('Results/Output_BUGS_Active_metrics_n_%d_p_%d_rep_%d_GuidanceKeep_%d.csv', n, p, rep, guidance_keep_val);
writetable(T_metrics, outfile_metrics);
fprintf('Saved compact metrics to: %s\n', outfile_metrics);

% ============================================================
% Save posterior beta draws for first true_p variables only
% ============================================================
beta_draws_top = out.beta_draws(:, 1:true_p);   % nSave x true_p

draw_id = (1:size(beta_draws_top,1))';
var_names = strcat("beta", string(1:true_p));

T_beta_draws = array2table(beta_draws_top, 'VariableNames', cellstr(var_names));
T_beta_draws = addvars(T_beta_draws, draw_id, 'Before', 1, 'NewVariableNames', 'draw_id');

outfile_draws = sprintf('Results/Output_BUGS_Active_beta_draws_top%d_n_%d_p_%d_rep_%d_GuidanceKeep_%d.csv', ...
    true_p, n, p, rep, guidance_keep_val);
writetable(T_beta_draws, outfile_draws);
fprintf('Saved posterior draws for first %d betas to: %s\n', true_p, outfile_draws);

% ============================================================
% Forest plot with posterior densities for first true_p variables
% ============================================================
true_p = 5;

beta_draws_top = out.beta_draws(:, 1:true_p);   % nSave x true_p
beta_mean_top  = mean(beta_draws_top, 1);
ci_2p5_top     = quantile(beta_draws_top, 0.025, 1);
ci_97p5_top    = quantile(beta_draws_top, 0.975, 1);

figure('Color', 'w', 'Position', [100, 100, 1050, 550]);
hold on;

ypos = true_p:-1:1;   % top-to-bottom ordering
half_width = 0.28;    % vertical thickness of density
x_all = beta_draws_top(:);
xpad = 0.08 * (max(x_all) - min(x_all) + eps);
xlim_vals = [min(x_all)-xpad, max(x_all)+ xpad];

for j = 1:true_p
    draws_j = beta_draws_top(:, j);

    % Kernel density
    [f, xi] = ksdensity(draws_j, 'NumPoints', 200);

    % Normalize density height for plotting
    f_scaled = f / max(f) * half_width;

    y0 = ypos(j);

    % Violin/density shape
    patch([xi, fliplr(xi)], ...
          [y0 + f_scaled, fliplr(y0 - f_scaled)], ...
          [0.75 0.82 0.92], ...
          'FaceAlpha', 0.85, ...
          'EdgeColor', [0.35 0.35 0.35], ...
          'LineWidth', 1.2);

    % Credible interval line
    plot([ci_2p5_top(j), ci_97p5_top(j)], [y0, y0], ...
        'k-', 'LineWidth', 3);

    beta_true_j = beta_true_std(j);   % make sure this is standardized

     plot(beta_true_j, y0, 'o', ...
    'MarkerSize', 10, ...
    'MarkerFaceColor', [0.85 0.33 0.10], ... % reddish fill
    'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5);
end

% Vertical reference line at 0
xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

% Axes / labels
xlim(xlim_vals);
ylim([0.5, true_p + 0.5]);

set(gca, ...
    'YTick', 1:true_p, ...
    'YTickLabel', compose('\\beta_{%d}', true_p:-1:1), ...
    'TickLabelInterpreter', 'tex', ...
    'FontSize', 18, ...
    'LineWidth', 1.2, ...
    'Box', 'on');
xlabel('Posterior coefficient', 'FontSize', 22);
ylabel('Signal', 'FontSize', 22);
title(sprintf('BUGS-Active (guidance: top %d)', guidance_keep_val), ...
    'FontSize', 25, 'FontWeight', 'bold');

grid on;
set(gca, 'GridAlpha', 0.2);

% ============================================================
% Metrics box (RHS outside plot)
% ============================================================

% Format text
txt = sprintf([ ...
    '\\bf Measures\\rm\n\n' ...
    'Time  : %.2f s\n' ...
    'RMSE(beta)  : %.3f\n' ...
    'MSE(y)  : %.3f\n' ...
    'TPR   : %.3f\n' ...
    'FPR   : %.3f\n' ...
    'FDR   : %.3f\n' ...
    'MCC   : %.3f\n'], ...
    metrics.CompTime, metrics.RMSE_beta, metrics.PredMSE, metrics.TPR, metrics.FPR, metrics.FDR, metrics.MCC);

% Shrink axes to leave space on RHS
ax = gca;
ax.Position = [0.10, 0.16, 0.65, 0.72];  
% [left, bottom, width, height]

% Add annotation box (normalized figure coordinates)
annotation('textbox', ...
    [0.76, 0.30, 0.20, 0.55], ...  % [x y width height]
    'String', txt, ...
    'Interpreter', 'tex', ...
    'FontSize', 16, ...
    'EdgeColor', 'black', ...
    'LineWidth', 1.4, ...
    'BackgroundColor', 'white', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle');


outfile_plot = sprintf('Results/Output_BUGS_Active_forest_density_top%d_n_%d_p_%d_rep_%d_GuidanceKeep_%d.png', ...
    true_p, n, p, rep, guidance_keep_val);
exportgraphics(gcf, outfile_plot, 'Resolution', 400);
fprintf('Saved forest-density plot to: %s\n', outfile_plot);