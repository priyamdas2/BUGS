clearvars;
addpath('./Data/');

% PARAMETERS
n = 500;
p = 1000000;
true_p = 10;
sigma = 1;
rho = 0.5;   
reps = 10;

% Ensure folder exists
outdir = 'Data';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% ---------------------------------------------------------
% Generate beta ONCE (fixed across reps)
% ---------------------------------------------------------
rng(1);  % fix beta reproducibly
[~, ~, beta] = generate_scenario_2_ultrahigh(n, true_p, p, sigma, rho, true);

% SAVE beta
fileB = sprintf('%s/Data_Scenario_2_beta_n_%d_p_%d.csv', outdir, n, p);
writematrix(beta, fileB);

% ---------------------------------------------------------
% Generate and save each replicate
% ---------------------------------------------------------
for r = 1:reps
    rng(r);  % reproducibility per replicate

    [X, y, ~] = generate_scenario_2_ultrahigh(n, true_p, p, sigma, rho, true);

    fileX = sprintf('%s/Data_Scenario_2_X_n_%d_p_%d_rep_%d.csv', outdir, n, p, r);
    fileY = sprintf('%s/Data_Scenario_2_Y_n_%d_p_%d_rep_%d.csv', outdir, n, p, r);

    writematrix(X, fileX);
    writematrix(y, fileY);
    
    fprintf('Saved replicate %d/%d\n', r, reps);
end

fprintf('Files saved in folder: %s\n', outdir);
fprintf('Beta file:\n%s\n', fileB);