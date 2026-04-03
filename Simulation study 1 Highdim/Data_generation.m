clearvars;
addpath('./Data/');

% PARAMETERS
n = 200;
p = 10000;
true_p = 10;
sigma = 1;
reps = 10;   % you said 10 files

% Ensure folder exists
outdir = 'Data';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Generate beta ONCE (fixed across reps)
rng(1);  % just to fix beta reproducibly
[~, ~, beta] = generate_scenario_1(n, true_p, p, sigma);

% SAVE beta
fileB = sprintf('%s/Data_Scenario_1_beta_n_%d_p_%d.csv', outdir, n, p);
writematrix(beta, fileB);

% Generate and save each replicate separately
for r = 1:reps
    rng(r);  % reproducibility per replicate

    [X, y, ~] = generate_scenario_1(n, true_p, p, sigma);

    fileX = sprintf('%s/Data_Scenario_1_X_n_%d_p_%d_rep_%d.csv', outdir, n, p, r);
    fileY = sprintf('%s/Data_Scenario_1_Y_n_%d_p_%d_rep_%d.csv', outdir, n, p, r);

    writematrix(X, fileX);
    writematrix(y, fileY);
end

fprintf('Files saved in folder: %s\n', outdir);
fprintf('Beta file:\n%s\n', fileB);