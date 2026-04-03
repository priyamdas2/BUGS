function GR_table = compute_GelmanRubin_BUGS(varargin)
% COMPUTEGELMANRUBIN_BUGS
% Compute Gelman-Rubin R-hat diagnostics from multiple BUGS chains.
%
% This function reads chain outputs saved by run_BUGS_multiple_chains_baseline.m
% and computes R-hat for selected hyperparameters and beta coordinates.
%
% USAGE:
%   GR_table = compute_GelmanRubin_BUGS();
%
%   GR_table = compute_GelmanRubin_BUGS( ...
%       'rootdir', 'U:\BUGS\Sensitivity and convergence', ...
%       'scenario', 1, ...
%       'rep', 1, ...
%       'n', 200, ...
%       'p', 1000, ...
%       'chain_ids', 1:4, ...
%       'output_csv', true);
%
% OUTPUT:
%   GR_table : MATLAB table with columns
%       Parameter
%       ParamType
%       Index
%       Rhat
%       NumChains
%       NumDrawsUsed
%
% NOTES:
%   - Requires at least 2 readable chain files.
%   - If chains have different numbers of saved draws, all are truncated to
%     the minimum common length before computing R-hat.
%   - Uses the classical Gelman-Rubin PSRF formula.

% ------------------------------------------------------------
% Parse inputs
% ------------------------------------------------------------
pr = inputParser;
addParameter(pr, 'rootdir', 'U:\BUGS\Sensitivity and convergence', @(x)ischar(x) || isstring(x));
addParameter(pr, 'scenario', 1, @(x)isnumeric(x) && isscalar(x));
addParameter(pr, 'rep', 1, @(x)isnumeric(x) && isscalar(x));
addParameter(pr, 'n', 200, @(x)isnumeric(x) && isscalar(x));
addParameter(pr, 'p', 1000, @(x)isnumeric(x) && isscalar(x));
addParameter(pr, 'chain_ids', 1:4, @(x)isnumeric(x) && isvector(x));
addParameter(pr, 'output_csv', true, @(x)islogical(x) && isscalar(x));
addParameter(pr, 'output_dir', '', @(x)ischar(x) || isstring(x));
addParameter(pr, 'output_file', '', @(x)ischar(x) || isstring(x));
parse(pr, varargin{:});

rootdir   = char(pr.Results.rootdir);
scenario  = pr.Results.scenario;
rep       = pr.Results.rep;
n         = pr.Results.n;
p         = pr.Results.p;
chain_ids = pr.Results.chain_ids;
output_csv = pr.Results.output_csv;
output_dir = char(pr.Results.output_dir);
output_file = char(pr.Results.output_file);

if isempty(output_dir)
    output_dir = fullfile(rootdir, 'Convergence');
end

if isempty(output_file)
    output_file = sprintf( ...
        'GelmanRubin_BUGS_scen_%d_rep_%d_n_%d_p_%d.csv', ...
        scenario, rep, n, p);
end

convdir = fullfile(rootdir, 'Convergence');

% ------------------------------------------------------------
% Fixed monitored indices
% ------------------------------------------------------------
signal_idx = [1, 5, 10];
noise_idx = unique([11, round((11 + p)/2), p]);

% ------------------------------------------------------------
% Read chain files
% ------------------------------------------------------------
chain_data = struct([]);
n_valid = 0;

fprintf('Reading BUGS chain files from:\n%s\n', convdir);

for cid = chain_ids
    matfile = fullfile(convdir, sprintf( ...
        'BUGS_chain_%d_scen_%d_rep_%d_n_%d_p_%d.mat', ...
        cid, scenario, rep, n, p));
    
    if ~isfile(matfile)
        warning('Missing chain file: %s', matfile);
        continue;
    end
    
    S = load(matfile);
    
    if ~isfield(S, 'out')
        warning('File does not contain variable ''out'': %s', matfile);
        continue;
    end
    
    out = S.out;
    
    required_fields = {'beta_draws', 'tau_draws', 'c2_draws', 'eta_draws', 'sigma2_draws'};
    missing_field = false;
    for k = 1:numel(required_fields)
        if ~isfield(out, required_fields{k})
            warning('Field %s missing in file %s', required_fields{k}, matfile);
            missing_field = true;
        end
    end
    if missing_field
        continue;
    end
    
    n_valid = n_valid + 1;
    chain_data(n_valid).chain_id = cid;
    chain_data(n_valid).file = matfile;
    chain_data(n_valid).beta_draws = out.beta_draws;
    chain_data(n_valid).tau_draws = out.tau_draws(:);
    chain_data(n_valid).c2_draws = out.c2_draws(:);
    chain_data(n_valid).eta_draws = out.eta_draws(:);
    chain_data(n_valid).sigma2_draws = out.sigma2_draws(:);
end

if n_valid < 2
    error('Need at least 2 valid chains to compute Gelman-Rubin diagnostics.');
end

fprintf('Using %d valid chains.\n', n_valid);

% ------------------------------------------------------------
% Determine common number of draws
% ------------------------------------------------------------
draw_counts = zeros(n_valid, 1);
for m = 1:n_valid
    draw_counts(m) = min([ ...
        size(chain_data(m).beta_draws, 1), ...
        length(chain_data(m).tau_draws), ...
        length(chain_data(m).c2_draws), ...
        length(chain_data(m).eta_draws), ...
        length(chain_data(m).sigma2_draws)]);
end

n_draws_common = min(draw_counts);

if n_draws_common < 2
    error('Common number of draws across chains is less than 2.');
end

fprintf('Using %d post-burnin kept draws per chain for R-hat.\n', n_draws_common);

% ------------------------------------------------------------
% Build summary rows
% ------------------------------------------------------------
rows = {};

% Hyperparameters
rows(end+1,:) = build_row('tau',    'hyperparameter', NaN, collect_scalar(chain_data, 'tau_draws', n_draws_common)); %#ok<AGROW>
rows(end+1,:) = build_row('c2',     'hyperparameter', NaN, collect_scalar(chain_data, 'c2_draws', n_draws_common)); %#ok<AGROW>
rows(end+1,:) = build_row('eta',    'hyperparameter', NaN, collect_scalar(chain_data, 'eta_draws', n_draws_common)); %#ok<AGROW>
rows(end+1,:) = build_row('sigma2', 'hyperparameter', NaN, collect_scalar(chain_data, 'sigma2_draws', n_draws_common)); %#ok<AGROW>

% Signal betas
for j = signal_idx
    pname = sprintf('beta_%d', j);
    rows(end+1,:) = build_row(pname, 'signal_beta', j, collect_beta(chain_data, j, n_draws_common)); %#ok<AGROW>
end

% Noise betas
for j = noise_idx
    pname = sprintf('beta_%d', j);
    rows(end+1,:) = build_row(pname, 'noise_beta', j, collect_beta(chain_data, j, n_draws_common)); %#ok<AGROW>
end

% ------------------------------------------------------------
% Convert to table
% ------------------------------------------------------------
n_rows = size(rows,1);
Parameter   = strings(n_rows,1);
ParamType   = strings(n_rows,1);
Index       = nan(n_rows,1);
Rhat        = nan(n_rows,1);
NumChains   = nan(n_rows,1);
NumDrawsUsed = nan(n_rows,1);

for i = 1:n_rows
    Parameter(i)    = string(rows{i,1});
    ParamType(i)    = string(rows{i,2});
    Index(i)        = rows{i,3};
    draws_mat       = rows{i,4};
    
    Rhat(i)         = gelman_rubin_rhat(draws_mat);
    NumChains(i)    = size(draws_mat, 2);
    NumDrawsUsed(i) = size(draws_mat, 1);
end

GR_table = table(Parameter, ParamType, Index, Rhat, NumChains, NumDrawsUsed);

% ------------------------------------------------------------
% Save CSV
% ------------------------------------------------------------
if output_csv
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    outcsv = fullfile(output_dir, output_file);
    writetable(GR_table, outcsv);
    fprintf('Gelman-Rubin table saved to:\n%s\n', outcsv);
end

% ------------------------------------------------------------
% Display
% ------------------------------------------------------------
disp(GR_table);

end

% ============================================================
% Helper functions
% ============================================================

function row = build_row(param_name, param_type, idx, draws_mat)
row = {param_name, param_type, idx, draws_mat};
end

function draws_mat = collect_scalar(chain_data, field_name, n_draws_common)
m = numel(chain_data);
draws_mat = nan(n_draws_common, m);
for k = 1:m
    x = chain_data(k).(field_name);
    draws_mat(:,k) = x(1:n_draws_common);
end
end

function draws_mat = collect_beta(chain_data, beta_index, n_draws_common)
m = numel(chain_data);
draws_mat = nan(n_draws_common, m);
for k = 1:m
    B = chain_data(k).beta_draws;
    if size(B,2) < beta_index
        error('beta_draws in chain %d has only %d columns; requested beta_%d.', ...
            chain_data(k).chain_id, size(B,2), beta_index);
    end
    draws_mat(:,k) = B(1:n_draws_common, beta_index);
end
end

function Rhat = gelman_rubin_rhat(draws_mat)
% draws_mat is n x m:
%   n = number of iterations per chain
%   m = number of chains

[n, m] = size(draws_mat);

if m < 2
    Rhat = NaN;
    return;
end

chain_means = mean(draws_mat, 1);
grand_mean  = mean(chain_means);

% Between-chain variance
B = n / (m - 1) * sum((chain_means - grand_mean).^2);

% Within-chain variance
chain_vars = var(draws_mat, 0, 1);  % sample variance within each chain
W = mean(chain_vars);

if W <= 0 || ~isfinite(W)
    Rhat = NaN;
    return;
end

% Marginal posterior variance estimate
Varhat = ((n - 1) / n) * W + (1 / n) * B;

if Varhat < 0 || ~isfinite(Varhat)
    Rhat = NaN;
    return;
end

Rhat = sqrt(Varhat / W);
end