clearvars;

% ============================================================
% Working directory / paths
% ============================================================

addpath('./Data/');
addpath('./supp funs/');

% ============================================================
% Scenario settings
% ============================================================

scenario = 1;
num_reps = 10;

n = 500;
p = 1000000;

% ============================================================
% Selection / tuning settings
% ============================================================

sel_thresh = 0.01;

% Number of CV folds for LASSO and SCAD
nfolds = 10;

% SCAD settings
scad_gamma = 3.7;
scad_num_lambda = 100;
scad_lambda_ratio = 1e-3;

scad_max_iter = 1000;
scad_tol = 1e-6;

% Create Results folder if needed
if ~exist('Results', 'dir')
    mkdir('Results');
end


% ============================================================
% Loop over simulation replicates
% ============================================================

parfor rep = 1:num_reps
    
    rng(rep);
    
    fprintf(['\nperforming data rep: %d | scenario: %d | ' ...
             'n: %d | p: %d\n'], ...
             rep, scenario, n, p);
    
    
    % ========================================================
    % File names
    % ========================================================
    
    fileX = sprintf( ...
        'Data_Scenario_%d_X_n_%d_p_%d_rep_%d.csv', ...
        scenario, n, p, rep);
    
    fileY = sprintf( ...
        'Data_Scenario_%d_Y_n_%d_p_%d_rep_%d.csv', ...
        scenario, n, p, rep);
    
    fileB = sprintf( ...
        'Data_Scenario_%d_beta_n_%d_p_%d.csv', ...
        scenario, n, p);
    
    
    % ========================================================
    % Read data
    % ========================================================
    
    X = readmatrix(fileX);
    y_raw = readmatrix(fileY);
    beta_true = readmatrix(fileB);
    
    % Force correct shapes
    y_raw = y_raw(:);
    beta_true = beta_true(:);
    
    
    % ========================================================
    % Basic checks
    % ========================================================
    
    if size(X,1) ~= n || size(X,2) ~= p
        error(['X has dimension %d x %d, expected ' ...
               '%d x %d.'], ...
               size(X,1), size(X,2), n, p);
    end
    
    if length(y_raw) ~= n
        error('y has length %d, expected %d.', ...
              length(y_raw), n);
    end
    
    if length(beta_true) ~= p
        error('beta_true has length %d, expected %d.', ...
              length(beta_true), p);
    end
    
    
    % ========================================================
    % Standardize y and corresponding beta_true
    % ========================================================
    
    y_raw = y_raw - mean(y_raw);
    sy = std(y_raw);
    
    if sy == 0
        error('Standard deviation of y_raw is zero.');
    end
    
    y = y_raw / sy;
    beta_true_std = beta_true / sy;
    
    
    % ========================================================
    % True-model fitted values
    %
    % beta_true is sparse, so avoid multiplying the entire
    % 500 x 1,000,000 matrix by beta_true.
    % ========================================================
    
    true_idx = find(beta_true_std ~= 0);
    
    pred_ref = X(:, true_idx) * beta_true_std(true_idx);
    
    
    % ========================================================
    % Step 1: Sure Independence Screening
    % ========================================================
    
    tic;
    
    % Fan and Lv (2008)
    d_sis = min(p, floor(n / log(n)));
    
    % Predictor standard deviations
    x_sd = std(X, 0, 1);
    
    % Absolute marginal correlation scores.
    %
    % y is already centered. Therefore X' y is unaffected by
    % centering each column of X. Division by sd(X_j) gives
    % the same ranking as absolute marginal correlation.
    sis_score = abs(X' * y) ./ x_sd';
    
    % Handle zero-variance / invalid predictors
    sis_score(~isfinite(sis_score)) = -Inf;
    
    % Rank predictors
    [~, sis_order] = sort(sis_score, 'descend');
    
    % Retain top d_sis
    sis_idx = sis_order(1:d_sis);
    
    % Reduced design matrix
    X_sis = X(:, sis_idx);
    
    sis_time = toc;
    
    fprintf(['  SIS retained %d of %d predictors | ' ...
             'time: %.3f sec\n'], ...
             d_sis, p, sis_time);
    
    
    % ========================================================
    % Step 2A: SIS-LASSO
    % ========================================================
    
    rng(rep);
    
    tic;
    
    [B_lasso, FitInfo_lasso] = lasso( ...
        X_sis, ...
        y, ...
        'Alpha', 1, ...
        'CV', nfolds, ...
        'Standardize', true);
    
    lasso_time = toc;
    
    
    % --------------------------------------------------------
    % Select lambda minimizing CV MSE
    % --------------------------------------------------------
    
    idx_lasso = FitInfo_lasso.IndexMinMSE;
    
    beta_lasso_reduced = B_lasso(:, idx_lasso);
    intercept_lasso = FitInfo_lasso.Intercept(idx_lasso);
    
    
    % --------------------------------------------------------
    % Expand back to original p-dimensional space
    % --------------------------------------------------------
    
    beta_hat_lasso = zeros(p, 1);
    beta_hat_lasso(sis_idx) = beta_lasso_reduced;
    
    
    % --------------------------------------------------------
    % Predictions
    %
    % No need to multiply X by the full million-dimensional
    % coefficient vector.
    % --------------------------------------------------------
    
    pred_est_lasso = ...
        intercept_lasso + X_sis * beta_lasso_reduced;
    
    
    % Total time = SIS + LASSO
    total_time_sis_lasso = sis_time + lasso_time;
    
    
    % --------------------------------------------------------
    % Metrics
    % --------------------------------------------------------
    
    metrics_lasso = compute_metrics_sis_method( ...
        beta_hat_lasso, ...
        beta_true_std, ...
        pred_est_lasso, ...
        pred_ref, ...
        total_time_sis_lasso, ...
        sel_thresh);
    
    
    % --------------------------------------------------------
    % Save SIS-LASSO
    % --------------------------------------------------------
    
    T_lasso = struct2table(metrics_lasso);
    
    outfile_lasso = sprintf( ...
        ['Results/Output_scen_%d_n_%d_p_%d_' ...
         'SIS_LASSO_rep_%d.csv'], ...
         scenario, n, p, rep);
    
    writetable(T_lasso, outfile_lasso);
    
    fprintf(['  SIS-LASSO fit time: %.3f sec | ' ...
             'total: %.3f sec\n'], ...
             lasso_time, total_time_sis_lasso);
    
    disp(T_lasso);
    
    fprintf('Saved metrics to: %s\n', outfile_lasso);
    
    
    % ========================================================
    % Step 2B: SIS-SCAD
    % ========================================================
    
    rng(rep);
    
    tic;
    
    [beta_scad_reduced, intercept_scad, lambda_scad] = ...
        fit_scad_cv( ...
            X_sis, ...
            y, ...
            nfolds, ...
            scad_gamma, ...
            rep, ...
            scad_num_lambda, ...
            scad_lambda_ratio, ...
            scad_max_iter, ...
            scad_tol);
    
    scad_time = toc;
    
    
    % --------------------------------------------------------
    % Expand back to original p-dimensional space
    % --------------------------------------------------------
    
    beta_hat_scad = zeros(p, 1);
    beta_hat_scad(sis_idx) = beta_scad_reduced;
    
    
    % --------------------------------------------------------
    % Predictions
    % --------------------------------------------------------
    
    pred_est_scad = ...
        intercept_scad + X_sis * beta_scad_reduced;
    
    
    % Total time = SIS + SCAD
    total_time_sis_scad = sis_time + scad_time;
    
    
    % --------------------------------------------------------
    % Metrics
    % --------------------------------------------------------
    
    metrics_scad = compute_metrics_sis_method( ...
        beta_hat_scad, ...
        beta_true_std, ...
        pred_est_scad, ...
        pred_ref, ...
        total_time_sis_scad, ...
        sel_thresh);
    
    
    % --------------------------------------------------------
    % Save SIS-SCAD
    % --------------------------------------------------------
    
    T_scad = struct2table(metrics_scad);
    
    outfile_scad = sprintf( ...
        ['Results/Output_scen_%d_n_%d_p_%d_' ...
         'SIS_SCAD_rep_%d.csv'], ...
         scenario, n, p, rep);
    
    writetable(T_scad, outfile_scad);
    
    fprintf(['  SIS-SCAD fit time: %.3f sec | ' ...
             'total: %.3f sec | lambda: %.6g\n'], ...
             scad_time, ...
             total_time_sis_scad, ...
             lambda_scad);
    
    disp(T_scad);
    
    fprintf('Saved metrics to: %s\n', outfile_scad);
    
end


% ========================================================================
% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================
% ========================================================================


function metrics = compute_metrics_sis_method( ...
    beta_est, beta_ref, pred_est, pred_ref, comp_time, tol)

% ============================================================
% Computes the same metrics used in the R SIS-LASSO/SCAD code.
%
% IMPORTANT:
% beta_est and beta_ref are both length p.
% Thus variables discarded by SIS are treated as estimated zero.
% ============================================================

beta_est = beta_est(:);
beta_ref = beta_ref(:);

% ------------------------------------------------------------
% Estimation metrics
% ------------------------------------------------------------

diff_beta = beta_est - beta_ref;

L2 = sqrt(sum(diff_beta.^2));
L1 = sum(abs(diff_beta));
RMSE = sqrt(mean(diff_beta.^2));

% ------------------------------------------------------------
% Prediction MSE
% ------------------------------------------------------------

PredMSE = mean((pred_ref - pred_est).^2);

% ------------------------------------------------------------
% Support recovery
% ------------------------------------------------------------

true_supp = (beta_ref ~= 0);
est_supp = (abs(beta_est) > tol);

TP = sum(true_supp & est_supp);
FP = sum(~true_supp & est_supp);
TN = sum(~true_supp & ~est_supp);
FN = sum(true_supp & ~est_supp);

TPR = TP / max(TP + FN, 1);
FPR = FP / max(FP + TN, 1);
FDR = FP / max(TP + FP, 1);

% ------------------------------------------------------------
% MCC
% ------------------------------------------------------------

mcc_num = double(TP) * double(TN) - ...
          double(FP) * double(FN);

aa = double(TP + FP);
bb = double(TP + FN);
cc = double(TN + FP);
dd = double(TN + FN);

mcc_den = sqrt(aa * bb * cc * dd);

if ~isfinite(mcc_den) || mcc_den <= 0
    MCC = NaN;
else
    MCC = mcc_num / mcc_den;
end

% ------------------------------------------------------------
% Output
% ------------------------------------------------------------

metrics.L2 = L2;
metrics.L1 = L1;
metrics.RMSE = RMSE;
metrics.PredMSE = PredMSE;

metrics.TPR = TPR;
metrics.FPR = FPR;
metrics.MCC = MCC;
metrics.FDR = FDR;

metrics.CompTime = comp_time;

metrics.TP = TP;
metrics.FP = FP;
metrics.TN = TN;
metrics.FN = FN;

metrics.sel_thresh = tol;

end


% ========================================================================
% SCAD WITH K-FOLD CROSS-VALIDATION
% ========================================================================

function [beta_best, intercept_best, lambda_best] = ...
    fit_scad_cv(X, y, K, gamma_scad, seed, ...
                num_lambda, lambda_ratio, max_iter, tol)

% ============================================================
% Gaussian SCAD regression with K-fold cross-validation.
%
% Objective:
%
%   (1/(2n)) ||y - X beta||^2
%       + sum_j p_lambda(|beta_j|)
%
% SCAD parameter gamma_scad is conventionally 3.7.
%
% A common lambda grid is used across folds.
% ============================================================

rng(seed);

[n, ~] = size(X);

y = y(:);


% ------------------------------------------------------------
% Construct lambda sequence
%
% Use standardized full data only to define the common grid.
% Each actual fold is standardized separately during fitting.
% ------------------------------------------------------------

muX = mean(X, 1);

X_centered = X - muX;

xscale = sqrt(sum(X_centered.^2, 1) / n);

xscale(~isfinite(xscale) | xscale <= 0) = 1;

X_std = X_centered ./ xscale;

y_centered = y - mean(y);

lambda_max = max(abs(X_std' * y_centered)) / n;

if ~isfinite(lambda_max) || lambda_max <= 0
    
    beta_best = zeros(size(X,2), 1);
    intercept_best = mean(y);
    lambda_best = 0;
    
    return;
end


lambda_min = lambda_max * lambda_ratio;

lambda_seq = exp( ...
    linspace( ...
        log(lambda_max), ...
        log(lambda_min), ...
        num_lambda));


% ------------------------------------------------------------
% Create deterministic balanced fold assignments
% ------------------------------------------------------------

perm = randperm(n);

fold_id = zeros(n, 1);

for ii = 1:n
    fold_id(perm(ii)) = mod(ii - 1, K) + 1;
end


% ------------------------------------------------------------
% Cross-validation
% ------------------------------------------------------------

cv_mse = zeros(K, num_lambda);

for k = 1:K
    
    train_idx = (fold_id ~= k);
    test_idx  = (fold_id == k);
    
    X_train = X(train_idx, :);
    y_train = y(train_idx);
    
    X_test = X(test_idx, :);
    y_test = y(test_idx);
    
    
    % Fit full SCAD path on training data
    [B_path, intercept_path] = scad_path_fit( ...
        X_train, ...
        y_train, ...
        lambda_seq, ...
        gamma_scad, ...
        max_iter, ...
        tol);
    
    
    % Predictions for every lambda
    pred = X_test * B_path;
    
    pred = pred + repmat( ...
        intercept_path, ...
        size(X_test,1), ...
        1);
    
    
    % Validation MSE for each lambda
    err = y_test - pred;
    
    cv_mse(k, :) = mean(err.^2, 1);
    
end


% ------------------------------------------------------------
% Lambda minimizing average CV error
% ------------------------------------------------------------

mean_cv_mse = mean(cv_mse, 1);

[~, best_idx] = min(mean_cv_mse);

lambda_best = lambda_seq(best_idx);


% ------------------------------------------------------------
% Refit full data across the common path and extract best lambda
% ------------------------------------------------------------

[B_full, intercept_full] = scad_path_fit( ...
    X, ...
    y, ...
    lambda_seq, ...
    gamma_scad, ...
    max_iter, ...
    tol);

beta_best = B_full(:, best_idx);
intercept_best = intercept_full(best_idx);

end


% ========================================================================
% FIT COMPLETE SCAD REGULARIZATION PATH
% ========================================================================

function [B_original, intercepts] = scad_path_fit( ...
    X, y, lambda_seq, gamma_scad, max_iter, tol)

[n, pp] = size(X);

y = y(:);

% ------------------------------------------------------------
% Standardize predictors
%
% Scale so that
%
%       (1/n) X_j' X_j = 1
%
% after centering.
% ------------------------------------------------------------

muX = mean(X, 1);

Xc = X - muX;

xscale = sqrt(sum(Xc.^2, 1) / n);

invalid_scale = ~isfinite(xscale) | xscale <= 0;
xscale(invalid_scale) = 1;

Z = Xc ./ xscale;

% Center response
muy = mean(y);
yc = y - muy;


% ------------------------------------------------------------
% Initialize
% ------------------------------------------------------------

num_lambda = length(lambda_seq);

B_standardized = zeros(pp, num_lambda);

beta = zeros(pp, 1);

resid = yc;


% ------------------------------------------------------------
% Follow lambda path from large -> small
%
% Warm start from previous lambda.
% ------------------------------------------------------------

for ll = 1:num_lambda
    
    lambda = lambda_seq(ll);
    
    for iter = 1:max_iter
        
        max_change = 0;
        
        for j = 1:pp
            
            beta_old = beta(j);
            
            % Add old contribution back into residual
            resid = resid + Z(:,j) * beta_old;
            
            % Because columns are standardized:
            %
            % (1/n) Z_j' Z_j = 1
            %
            z_j = (Z(:,j)' * resid) / n;
            
            % SCAD coordinate update
            beta_new = scad_threshold( ...
                z_j, lambda, gamma_scad);
            
            beta(j) = beta_new;
            
            % Remove new contribution
            resid = resid - Z(:,j) * beta_new;
            
            change = abs(beta_new - beta_old);
            
            if change > max_change
                max_change = change;
            end
            
        end
        
        if max_change < tol
            break;
        end
        
    end
    
    B_standardized(:, ll) = beta;
    
end


% ------------------------------------------------------------
% Transform coefficients back to original predictor scale
% ------------------------------------------------------------

B_original = B_standardized ./ xscale(:);

% Zero-variance predictors must remain zero
B_original(invalid_scale, :) = 0;

% Corresponding intercepts
intercepts = muy - muX * B_original;

end


% ========================================================================
% SCAD COORDINATE THRESHOLDING OPERATOR
% ========================================================================

function beta = scad_threshold(z, lambda, a)

az = abs(z);

if az <= lambda
    
    beta = 0;
    
elseif az <= 2 * lambda
    
    beta = sign(z) * (az - lambda);
    
elseif az <= a * lambda
    
    beta = ((a - 1) * z - ...
            sign(z) * a * lambda) / ...
            (a - 2);
    
else
    
    beta = z;
    
end

end