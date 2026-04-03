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

% ============================================================
% Fixed convergence-analysis setup
% ============================================================
scenario = 1;
rep = 1;                 % use only first replicate
n = 200;
p = 1000;
true_p = 10;

nchains = 4;

outdir = fullfile(rootdir, 'Convergence');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

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
% Standardize y exactly as in simulation code
% ============================================================
y_raw = y_raw - mean(y_raw);
sy = std(y_raw);

if sy == 0
    error('Standard deviation of y_raw is zero.');
end

y = y_raw / sy;
beta_true_std = beta_true / sy;

% ============================================================
% Baseline BUGS options
% ============================================================
base_opts = struct();

base_opts.niter = 5000;
base_opts.burnin = 1000;
base_opts.thin = 5;
base_opts.display_every = 100;
base_opts.verbose = true;

base_opts.tau0 = 0.001;
base_opts.sigma_eta = 0.002;
base_opts.Cz = 2;
base_opts.a_c = 10;
base_opts.b_c = 100;
base_opts.a_sigma = 2;
base_opts.b_sigma = 2;

% Optional sampler controls
base_opts.eps_score = 1e-6;
base_opts.slice_w = 0.5;
base_opts.slice_m = 100;
base_opts.save_lambda_max = 100;

% ============================================================
% Fixed monitored indices for later convergence summaries
% ============================================================
signal_idx = [1, 5, 10];
noise_idx = unique([11, round((11 + p)/2), p]);

fprintf('Signal indices monitored: %s\n', mat2str(signal_idx));
fprintf('Noise indices monitored : %s\n', mat2str(noise_idx));

% ============================================================
% Crude data-informed beta start for one chain
% ============================================================
beta_marg = (X' * y) / n;
beta_marg = beta_marg(:);

% ============================================================
% Run chains
% ============================================================
for chain_id = 1:nchains
    
    opts = base_opts;
    opts.seed = 100 + chain_id;   % distinct RNG seeds
    
    % --------------------------------------------------------
    % Dispersed initial values across chains
    % --------------------------------------------------------
    init = struct();
    
    switch chain_id
        
        case 1
            % Conservative / tight shrinkage start
            init.beta = zeros(p,1);
            init.lambda2 = ones(p,1);
            init.tau2 = opts.tau0^2 / 100;
            init.c2 = 1;
            init.eta = 0.01;
            init.sigma2 = 1;
            init_label = 'conservative';
            
        case 2
            % Default-ish start
            init.beta = zeros(p,1);
            init.lambda2 = ones(p,1);
            init.tau2 = opts.tau0^2 / 10;
            init.c2 = 1;
            init.eta = 0.1;
            init.sigma2 = 1;
            init_label = 'defaultish';
            
        case 3
            % More permissive start
            init.beta = zeros(p,1);
            init.lambda2 = 5 * ones(p,1);
            init.tau2 = opts.tau0^2;
            init.c2 = 10;
            init.eta = 0.5;
            init.sigma2 = 2;
            init_label = 'permissive';
            
        case 4
            % Crude data-informed start
            init.beta = 0.1 * beta_marg;
            init.lambda2 = ones(p,1);
            init.tau2 = 5 * opts.tau0^2;
            init.c2 = 50;
            init.eta = 0.2;
            init.sigma2 = 1;
            init_label = 'data_informed';
            
        otherwise
            error('Unsupported chain_id = %d', chain_id);
    end
    
    opts.init = init;
    
    fprintf('\n============================================\n');
    fprintf('Running chain %d / %d (%s)\n', chain_id, nchains, init_label);
    fprintf('Seed = %d\n', opts.seed);
    fprintf('Initial tau = %.6f\n', sqrt(init.tau2));
    fprintf('Initial c2  = %.6f\n', init.c2);
    fprintf('Initial eta = %.6f\n', init.eta);
    fprintf('Initial sigma2 = %.6f\n', init.sigma2);
    fprintf('============================================\n');
    
    % --------------------------------------------------------
    % Run BUGS chain
    % --------------------------------------------------------
    out = BUGS_mcmc_init(X, y, opts);
    
    % --------------------------------------------------------
    % Save chain output
    % --------------------------------------------------------
    savefile = fullfile(outdir, sprintf( ...
        'BUGS_chain_%d_scen_%d_rep_%d_n_%d_p_%d.mat', ...
        chain_id, scenario, rep, n, p));
    
    chain_info = struct();
    chain_info.chain_id = chain_id;
    chain_info.init_label = init_label;
    chain_info.seed = opts.seed;
    chain_info.scenario = scenario;
    chain_info.rep = rep;
    chain_info.n = n;
    chain_info.p = p;
    chain_info.true_p = true_p;
    chain_info.signal_idx = signal_idx;
    chain_info.noise_idx = noise_idx;
    
    save(savefile, 'out', 'chain_info', 'beta_true', 'beta_true_std', '-v7.3');
    
    fprintf('Saved chain %d output to:\n%s\n', chain_id, savefile);
end

fprintf('\nAll %d BUGS chains completed successfully.\n', nchains);