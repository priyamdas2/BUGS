% ============================================================
% Prepare raw-data files for final prediction plot
% Reads:
%   - full raw X from X_matrix_final.mat
%   - full raw y from Pheno_for_MATLAB.csv
%   - top-ranked indices from Output_FULL_topK_beta_summary.csv
%
% Saves into folder:
%   ./Final prediction plot data/
%
% Outputs:
%   rawX_mean.csv
%   rawX_sd.csv
%   rawY_mean_sd.csv
%   selected_X_raw_topK.csv
%   y_raw.csv
%   selected_topK_metadata.csv
%   prediction_plot_data.mat
% ============================================================

clearvars -except X_raw_FULL y_raw_FULL
clc;

addpath('./BUGS/');
addpath('./Dummy data/');
addpath('./Age methylation data/');

%% =========================
% User settings
% ==========================
fileX       = 'X_matrix_final.mat';
fileY       = 'Pheno_for_MATLAB.csv';
fileSummary = './Output/Output_FULL_topK_beta_summary.csv';

outdir = './Final prediction plot data/';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% choose how many top-ranked CpGs to extract
topK_extract = 200;   % e.g. 10, 20, 50, 200

% response column name in phenotype CSV
y_var_name = 'Age_months';

%% =========================
% Read full raw X
% ==========================
if ~exist('X_raw_FULL','var') || isempty(X_raw_FULL)
    fprintf('Loading full raw X...\n');
    t_loadX = tic;
    S = load(fileX);
    fn = fieldnames(S);
    X_raw_FULL = S.(fn{1});
    clear S fn;
    fprintf('Time to load full X: %.2f seconds\n', toc(t_loadX));
else
    fprintf('Full raw X already loaded. Skipping load.\n');
end

%% =========================
% Read full raw y
% ==========================
if ~exist('y_raw_FULL','var') || isempty(y_raw_FULL)
    fprintf('Loading full raw y...\n');
    Ty = readtable(fileY);

    if ~ismember(y_var_name, Ty.Properties.VariableNames)
        error('Column "%s" not found in %s.', y_var_name, fileY);
    end

    y_raw_FULL = Ty.(y_var_name);
    y_raw_FULL = y_raw_FULL(:);
    clear Ty;
else
    fprintf('Full raw y already loaded. Skipping load.\n');
end

%% =========================
% Basic checks
% ==========================
[n_full, p_full] = size(X_raw_FULL);

fprintf('Loaded raw data with n = %d and p = %d.\n', n_full, p_full);

if length(y_raw_FULL) ~= n_full
    error('y has length %d, but X has %d rows.', length(y_raw_FULL), n_full);
end

%% =========================
% Compute raw-data summaries
% ==========================
fprintf('Computing raw X / y summaries...\n');

rawX_mean = mean(X_raw_FULL, 1);
rawX_sd   = std(X_raw_FULL, 0, 1);

rawY_mean = mean(y_raw_FULL);
rawY_sd   = std(y_raw_FULL);

% protect against zero-variance columns if needed downstream
rawX_sd_safe = rawX_sd;
rawX_sd_safe(rawX_sd_safe == 0) = 1;

%% =========================
% Read top-ranked summary file
% ==========================
fprintf('Reading top-ranked summary file...\n');

Tsum = readtable(fileSummary);

reqVars = {'rank', 'index'};
for j = 1:length(reqVars)
    if ~ismember(reqVars{j}, Tsum.Properties.VariableNames)
        error('Required column "%s" not found in %s.', reqVars{j}, fileSummary);
    end
end

% sort by rank and keep top K
Tsum = sortrows(Tsum, 'rank', 'ascend');
topK_extract = min(topK_extract, height(Tsum));
Ttop = Tsum(1:topK_extract, :);

selected_idx = Ttop.index(:);

% sanity checks
if any(selected_idx < 1) || any(selected_idx > p_full)
    error('Some selected indices are outside valid column range 1..%d.', p_full);
end

if numel(unique(selected_idx)) < numel(selected_idx)
    warning('Duplicate predictor indices found among selected top ranks.');
end

fprintf('Extracting top %d raw X columns based on rank.\n', topK_extract);

%% =========================
% Extract selected raw X
% ==========================
X_raw_selected = X_raw_FULL(:, selected_idx);

%% =========================
% Build metadata table
% ==========================
selected_rank = Ttop.rank(:);
selected_index = Ttop.index(:);

Meta = table(selected_rank, selected_index, ...
    'VariableNames', {'rank', 'index'});

% add optional summary columns if present
optionalVars = {'beta_mean', 'beta_median', 'beta_sd', 'ci_2p5', 'ci_97p5', ...
                'beta_mean_raw', 'beta_median_raw', 'ci_2p5_raw', 'ci_97p5_raw', ...
                'post_prob_abs_gt_thresh', 'selected_mean_thresh', ...
                'selected_median_thresh', 'selected_ci_excludes_zero', ...
                'sign_beta_mean', 'sign_beta_median'};

for j = 1:length(optionalVars)
    v = optionalVars{j};
    if ismember(v, Ttop.Properties.VariableNames)
        Meta.(v) = Ttop.(v);
    end
end

% also store corresponding raw X mean and sd for the selected columns
Meta.rawX_mean = rawX_mean(selected_idx).';
Meta.rawX_sd   = rawX_sd(selected_idx).';

%% =========================
% Write outputs
% ==========================
fprintf('Writing outputs to %s\n', outdir);

% full raw X summaries
writematrix(rawX_mean(:), fullfile(outdir, 'rawX_mean.csv'));
writematrix(rawX_sd(:),   fullfile(outdir, 'rawX_sd.csv'));

% raw y summaries
TySummary = table(rawY_mean, rawY_sd, ...
    'VariableNames', {'rawY_mean', 'rawY_sd'});
writetable(TySummary, fullfile(outdir, 'rawY_mean_sd.csv'));

% selected raw X and full raw y
writematrix(X_raw_selected, fullfile(outdir, sprintf('selected_X_raw_top%d.csv', topK_extract)));
writematrix(y_raw_FULL(:),  fullfile(outdir, 'y_raw.csv'));

% metadata for selected columns
writetable(Meta, fullfile(outdir, sprintf('selected_top%d_metadata.csv', topK_extract)));

% save all in one MAT file too
save(fullfile(outdir, sprintf('prediction_plot_data_top%d.mat', topK_extract)), ...
    'X_raw_selected', 'y_raw_FULL', ...
    'rawX_mean', 'rawX_sd', 'rawX_sd_safe', ...
    'rawY_mean', 'rawY_sd', ...
    'selected_idx', 'Meta', 'topK_extract', '-v7.3');

%% =========================
% Console summary
% ==========================
fprintf('\nDone.\n');
fprintf('Saved files:\n');
fprintf('  %s\n', fullfile(outdir, 'rawX_mean.csv'));
fprintf('  %s\n', fullfile(outdir, 'rawX_sd.csv'));
fprintf('  %s\n', fullfile(outdir, 'rawY_mean_sd.csv'));
fprintf('  %s\n', fullfile(outdir, sprintf('selected_X_raw_top%d.csv', topK_extract)));
fprintf('  %s\n', fullfile(outdir, 'y_raw.csv'));
fprintf('  %s\n', fullfile(outdir, sprintf('selected_top%d_metadata.csv', topK_extract)));
fprintf('  %s\n', fullfile(outdir, sprintf('prediction_plot_data_top%d.mat', topK_extract)));

fprintf('\nFirst few selected indices by rank:\n');
disp(Meta(1:min(10,height(Meta)), :));