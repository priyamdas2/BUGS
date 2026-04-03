clearvars;
addpath('./Data/');

% PARAMETERS
n = 200;
p = 100000;   % change to 100000 for ultra-high
true_p = 10;
sigma = 1;
reps = 10;   

% Ultra-high settings
use_single = true;     % use single precision for speed + memory
block_size = 5000;     % tune if needed (e.g., 2000–10000)

% Ensure folder exists
outdir = 'Data';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Generate beta ONCE (fixed across reps)
rng(1);  % just to fix beta reproducibly
[~, ~, beta] = generate_scenario_1_ultrahigh(n, true_p, p, sigma, use_single, block_size);

% SAVE beta
fileB = sprintf('%s/Data_Scenario_1_beta_n_%d_p_%d.csv', outdir, n, p);
writematrix(beta, fileB);

% Generate and save each replicate separately
for r = 1:reps
    rng(r);  % reproducibility per replicate

    [X, y, ~] = generate_scenario_1_ultrahigh(n, true_p, p, sigma, use_single, block_size);

    fileX = sprintf('%s/Data_Scenario_1_X_n_%d_p_%d_rep_%d.csv', outdir, n, p, r);
    fileY = sprintf('%s/Data_Scenario_1_Y_n_%d_p_%d_rep_%d.csv', outdir, n, p, r);

    % Write to CSV (note: still I/O heavy for large p)
    writematrix(X, fileX);
    writematrix(y, fileY);

    fprintf('Saved replicate %d/%d\n', r, reps);
end

fprintf('Files saved in folder: %s\n', outdir);
fprintf('Beta file:\n%s\n', fileB);