clearvars;
clc;
close all;

% ============================================================
% Paths
% ============================================================
rootdir = 'U:\BUGS\Sensitivity and convergence';
cd(rootdir);

addpath(fullfile(rootdir, 'BUGS'));
addpath(fullfile(rootdir, 'Data'));
addpath(fullfile(rootdir, 'supp funs'));

results_dir = fullfile(rootdir, 'Results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% ============================================================
% Fixed setup: one default BUGS run
% ============================================================
scenario = 1;
rep = 1;
n = 200;
p = 1000;

% ============================================================
% File names
% ============================================================
fileX = fullfile(rootdir, 'Data', ...
    sprintf('Data_Scenario_%d_X_n_%d_p_%d_rep_%d.csv', scenario, n, p, rep));

fileY = fullfile(rootdir, 'Data', ...
    sprintf('Data_Scenario_%d_Y_n_%d_p_%d_rep_%d.csv', scenario, n, p, rep));

fileB = fullfile(rootdir, 'Data', ...
    sprintf('Data_Scenario_%d_beta_n_%d_p_%d.csv', scenario, n, p));

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
% Standardize y exactly as in your simulation code
% ============================================================
y_raw = y_raw - mean(y_raw);
sy = std(y_raw);

if sy == 0
    error('Standard deviation of y_raw is zero.');
end

y = y_raw / sy;
beta_true_std = beta_true / sy;

% ============================================================
% Default BUGS options
% ============================================================
opts = struct();
opts.niter = 5000;
opts.burnin = 1000;
opts.thin = 5;
opts.display_every = 100;
opts.verbose = true;

opts.tau0 = 0.001;
opts.sigma_eta = 0.002;
opts.Cz = 2;
opts.a_c = 10;
opts.b_c = 100;
opts.a_sigma = 2;
opts.b_sigma = 2;

opts.eps_score = 1e-6;
opts.slice_w = 0.5;
opts.slice_m = 100;
opts.save_lambda_max = 100;
opts.seed = 1;

fprintf('Running default BUGS diagnostics fit: Scenario %d, Rep %d, n=%d, p=%d\n', ...
    scenario, rep, n, p);
disp(opts);

% ============================================================
% Fit BUGS once
% ============================================================
tic;
out = BUGS_mcmc(X, y, opts);
comp_time = toc;

% Save raw fit
fit_matfile = fullfile(results_dir, sprintf( ...
    'BUGS_default_fit_scen_%d_rep_%d_n_%d_p_%d.mat', ...
    scenario, rep, n, p));
save(fit_matfile, 'out', 'opts', 'beta_true', 'beta_true_std', 'comp_time', '-v7.3');

% ============================================================
% Extract saved draws
% ============================================================
beta_draws = out.beta_draws;
tau_draws = out.tau_draws(:);
c2_draws = out.c2_draws(:);
eta_draws = out.eta_draws(:);
sigma2_draws = out.sigma2_draws(:);

nkeep = size(beta_draws, 1);

% ============================================================
% Monitored parameters
% ============================================================
signal_idx = [1, 5, 10];
noise_idx = [11, 506, 1000];

param_names = { ...
    'tau', 'c2', ...
    'eta', 'sigma2', ...
    'beta_1', 'beta_5', ...
    'beta_10', 'beta_11', ...
    'beta_506', 'beta_1000'};

param_types = { ...
    'hyperparameter', 'hyperparameter', ...
    'hyperparameter', 'hyperparameter', ...
    'signal_beta', 'signal_beta', ...
    'signal_beta', 'noise_beta', ...
    'noise_beta', 'noise_beta'};

param_index = [ ...
    NaN, NaN, ...
    NaN, NaN, ...
    1, 5, ...
    10, 11, ...
    506, 1000];

draw_list = cell(10,1);
draw_list{1} = tau_draws;
draw_list{2} = c2_draws;
draw_list{3} = eta_draws;
draw_list{4} = sigma2_draws;
draw_list{5} = beta_draws(:,1);
draw_list{6} = beta_draws(:,5);
draw_list{7} = beta_draws(:,10);
draw_list{8} = beta_draws(:,11);
draw_list{9} = beta_draws(:,506);
draw_list{10} = beta_draws(:,1000);

% ============================================================
% Figure style
% ============================================================
fig_w = 1400;
fig_h = 1500;

lw_trace = 1.2;
lw_mean  = 1.2;
lw_acf   = 1.0;
lw_ref   = 0.8;

fs_label = 13;    % axis labels
fs_sg    = 18;    % super title
fs_axis  = 13;    % tick labels
fs_title = 16;    % subplot titles

% LaTeX-style parameter labels
param_labels = { ...
    '$\tau$', '$c^2$', ...
    '$\eta$', '$\sigma^2$', ...
    '$\beta_1$', '$\beta_5$', ...
    '$\beta_{10}$', '$\beta_{11}$', ...
    '$\beta_{506}$', '$\beta_{1000}$'};

% ============================================================
% 1. Trace plots
% ============================================================
f1 = figure('Color', 'w', 'Position', [50, 50, fig_w, fig_h]);

for i = 1:10
    subplot(5,2,i);
    x = draw_list{i};

    plot(1:nkeep, x, 'k-', 'LineWidth', lw_trace);
    box off;
    grid on;

    ax = gca;
    ax.FontSize = fs_axis;
    ax.LineWidth = 1;

    title(param_labels{i}, ...
        'Interpreter', 'latex', ...
        'FontSize', fs_title, ...
        'FontWeight', 'bold');

    xlabel('Iteration', 'FontSize', fs_label);
    ylabel('');
end

sgtitle(sprintf('BUGS trace plots ($n = %d$, $p = %d$)', n, p), ...
    'Interpreter', 'latex', ...
    'FontSize', fs_sg, ...
    'FontWeight', 'bold');

trace_file_png = fullfile(results_dir, sprintf( ...
    'BUGS_traceplots_default_scen_%d_rep_%d_n_%d_p_%d.png', ...
    scenario, rep, n, p));
trace_file_pdf = fullfile(results_dir, sprintf( ...
    'BUGS_traceplots_default_scen_%d_rep_%d_n_%d_p_%d.pdf', ...
    scenario, rep, n, p));

exportgraphics(f1, trace_file_png, 'Resolution', 400);
exportgraphics(f1, trace_file_pdf, 'ContentType', 'vector');

% ============================================================
% 2. Running means
% ============================================================
f2 = figure('Color', 'w', 'Position', [100, 60, fig_w, fig_h]);

for i = 1:10
    subplot(5,2,i);
    x = draw_list{i};
    running_mean = cumsum(x) ./ (1:numel(x))';

    plot(1:nkeep, running_mean, 'k-', 'LineWidth', lw_mean);
    box off;
    grid on;

    ax = gca;
    ax.FontSize = fs_axis;
    ax.LineWidth = 1;

    title(param_labels{i}, ...
        'Interpreter', 'latex', ...
        'FontSize', fs_title, ...
        'FontWeight', 'bold');

    xlabel('Iteration', 'FontSize', fs_label);
    ylabel('');
end

sgtitle(sprintf('BUGS running means ($n = %d$, $p = %d$)', n, p), ...
    'Interpreter', 'latex', ...
    'FontSize', fs_sg, ...
    'FontWeight', 'bold');

runmean_file_png = fullfile(results_dir, sprintf( ...
    'BUGS_runningmeans_default_scen_%d_rep_%d_n_%d_p_%d.png', ...
    scenario, rep, n, p));
runmean_file_pdf = fullfile(results_dir, sprintf( ...
    'BUGS_runningmeans_default_scen_%d_rep_%d_n_%d_p_%d.pdf', ...
    scenario, rep, n, p));

exportgraphics(f2, runmean_file_png, 'Resolution', 400);
exportgraphics(f2, runmean_file_pdf, 'ContentType', 'vector');

% ============================================================
% 3. Autocorrelation plots
% ============================================================
max_lag = min(50, floor(nkeep/4));

f3 = figure('Color', 'w', 'Position', [150, 70, fig_w, fig_h]);

for i = 1:10
    subplot(5,2,i);
    x = draw_list{i};
    acf_vals = compute_acf(x, max_lag);

    stem(0:max_lag, acf_vals, 'k', 'LineWidth', lw_acf, 'Marker', 'none');
    hold on;
    conf = 1.96 / sqrt(numel(x));
    yline(conf, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', lw_ref);
    yline(-conf, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', lw_ref);
    yline(0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', lw_ref);
    hold off;

    box off;
    grid on;

    ax = gca;
    ax.FontSize = fs_axis;
    ax.LineWidth = 1;
    xlim([0 max_lag]);

    title(param_labels{i}, ...
        'Interpreter', 'latex', ...
        'FontSize', fs_title, ...
        'FontWeight', 'bold');

    xlabel('Lag', 'FontSize', fs_label);
    ylabel('');
end

sgtitle(sprintf('BUGS autocorrelation plots ($n = %d$, $p = %d$)', n, p), ...
    'Interpreter', 'latex', ...
    'FontSize', fs_sg, ...
    'FontWeight', 'bold');

acf_file_png = fullfile(results_dir, sprintf( ...
    'BUGS_acf_default_scen_%d_rep_%d_n_%d_p_%d.png', ...
    scenario, rep, n, p));
acf_file_pdf = fullfile(results_dir, sprintf( ...
    'BUGS_acf_default_scen_%d_rep_%d_n_%d_p_%d.pdf', ...
    scenario, rep, n, p));

exportgraphics(f3, acf_file_png, 'Resolution', 400);
exportgraphics(f3, acf_file_pdf, 'ContentType', 'vector');

% ============================================================
% 4. Posterior density plots
% ============================================================
f4 = figure('Color', 'w', 'Position', [200, 80, fig_w, fig_h]);

for i = 1:10
    subplot(5,2,i);
    x = draw_list{i};
    [f, xi] = ksdensity(x);

    plot(xi, f, 'k-', 'LineWidth', lw_mean);
    box off;
    grid on;

    ax = gca;
    ax.FontSize = fs_axis;
    ax.LineWidth = 1;

    title(param_labels{i}, ...
        'Interpreter', 'latex', ...
        'FontSize', fs_title, ...
        'FontWeight', 'bold');

    xlabel('Value', 'FontSize', fs_label);
    ylabel('');
end

sgtitle(sprintf('BUGS posterior densities ($n = %d$, $p = %d$)', n, p), ...
    'Interpreter', 'latex', ...
    'FontSize', fs_sg, ...
    'FontWeight', 'bold');

dens_file_png = fullfile(results_dir, sprintf( ...
    'BUGS_density_default_scen_%d_rep_%d_n_%d_p_%d.png', ...
    scenario, rep, n, p));
dens_file_pdf = fullfile(results_dir, sprintf( ...
    'BUGS_density_default_scen_%d_rep_%d_n_%d_p_%d.pdf', ...
    scenario, rep, n, p));

exportgraphics(f4, dens_file_png, 'Resolution', 400);
exportgraphics(f4, dens_file_pdf, 'ContentType', 'vector');

% ============================================================
% 5. ESS and posterior summary table
% ============================================================
ESS = nan(10,1);
PostMean = nan(10,1);
PostMedian = nan(10,1);
PostSD = nan(10,1);
CI2p5 = nan(10,1);
CI97p5 = nan(10,1);

for i = 1:10
    x = draw_list{i};
    ESS(i) = compute_ess(x);
    PostMean(i) = mean(x);
    PostMedian(i) = median(x);
    PostSD(i) = std(x);
    CI2p5(i) = quantile(x, 0.025);
    CI97p5(i) = quantile(x, 0.975);
end

DiagnosticsTable = table( ...
    string(param_names(:)), ...
    string(param_types(:)), ...
    param_index(:), ...
    ESS, PostMean, PostMedian, PostSD, CI2p5, CI97p5, ...
    repmat(rep, 10, 1), ...
    repmat(nkeep, 10, 1), ...
    repmat(comp_time, 10, 1), ...
    'VariableNames', { ...
    'Parameter', 'ParamType', 'Index', ...
    'ESS', 'PostMean', 'PostMedian', 'PostSD', 'CI2p5', 'CI97p5', ...
    'Rep', 'NumDrawsUsed', 'CompTime'});

diag_csv = fullfile(results_dir, sprintf( ...
    'BUGS_singlechain_diagnostics_table_default_scen_%d_rep_%d_n_%d_p_%d.csv', ...
    scenario, rep, n, p));

writetable(DiagnosticsTable, diag_csv);

disp(DiagnosticsTable);

fprintf('\nSaved files:\n');
fprintf('%s\n', fit_matfile);
fprintf('%s\n', trace_file_png);
fprintf('%s\n', trace_file_pdf);
fprintf('%s\n', runmean_file_png);
fprintf('%s\n', runmean_file_pdf);
fprintf('%s\n', acf_file_png);
fprintf('%s\n', acf_file_pdf);
fprintf('%s\n', dens_file_png);
fprintf('%s\n', dens_file_pdf);
fprintf('%s\n', diag_csv);

% ============================================================
% Local helper functions
% ============================================================
function acf_vals = compute_acf(x, max_lag)
    x = x(:);
    x = x - mean(x);
    n = length(x);
    v = sum(x.^2) / n;
    acf_vals = zeros(max_lag+1, 1);
    for lag = 0:max_lag
        acf_vals(lag+1) = sum(x(1:n-lag) .* x(1+lag:n)) / n / v;
    end
end

function ess = compute_ess(x)
    x = x(:);
    n = length(x);
    max_lag = min(n-1, floor(n/2));
    acf_vals = compute_acf(x, max_lag);

    s = 0;
    for lag = 1:max_lag
        if acf_vals(lag+1) <= 0
            break;
        end
        s = s + acf_vals(lag+1);
    end

    tau_int = 1 + 2*s;
    ess = n / tau_int;
    ess = min(ess, n);
end