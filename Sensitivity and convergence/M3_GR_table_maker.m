clearvars;
clc;

% ============================================================
% Paths
% ============================================================
rootdir = 'U:\BUGS\Sensitivity and convergence';
cd(rootdir);

addpath(fullfile(rootdir, 'BUGS'));
addpath(fullfile(rootdir, 'Data'));
addpath(fullfile(rootdir, 'supp funs'));

GR_table = compute_GelmanRubin_BUGS();

% ============================================================
% Save to Results folder
% ============================================================
results_dir = fullfile(rootdir, 'Results');

if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

outfile = fullfile(results_dir, 'ALL_GR_table.csv');
writetable(GR_table, outfile);

fprintf('Saved Gelman-Rubin table to:\n%s\n', outfile);