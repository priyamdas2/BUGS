% ============================================================
% Simulation study with BayesReg (Lasso / Horseshoe / Horseshoe+)
% ============================================================
clearvars; clc; close all;
addpath('./supp funs/');
addpath('./Data/');
addpath('./bayesreg/');

scenario = 2;
num_reps = 10;

n = 200;
p = 1000;

burn = 1000;
iter = 5000;


for rep = 1:num_reps

    fprintf('performing data rep: %d | scenario: %d | n: %d | p: %d\n', rep, scenario, n, p);
    
    % -----------------------------
    % Toggle which methods to run
    % -----------------------------
    run_bayeslasso     = 1;
    run_horseshoe      = 1;
    run_horseshoeplus  = 1;
    
    if (run_bayeslasso + run_horseshoe + run_horseshoeplus) == 0
        error('Set at least one method to 1.');
    end
    
    % -----------------------------
    % User settings
    % -----------------------------
    
    thin = 1;
    
    rng(rep);
    
    use_topk = false;
    tol = 0.01;
    
    data_dir = 'Data';
    results_dir = 'Results';
    
    if ~exist(results_dir, 'dir')
        mkdir(results_dir);
    end
    
    % -----------------------------
    % Read data
    % -----------------------------
    fileX = sprintf('%s/Data_Scenario_%d_X_n_%d_p_%d_rep_%d.csv', ...
        data_dir, scenario, n, p, rep);
    fileY = sprintf('%s/Data_Scenario_%d_Y_n_%d_p_%d_rep_%d.csv', ...
        data_dir, scenario, n, p, rep);
    fileB = sprintf('%s/Data_Scenario_%d_beta_n_%d_p_%d.csv', ...
        data_dir, scenario, n, p);
    
    X = readmatrix(fileX);
    y_raw = readmatrix(fileY);
    beta_true = readmatrix(fileB);
    
    y_raw = y_raw(:);
    beta_true = beta_true(:);
    
    % -----------------------------
    % Standardize y and beta_true
    % -----------------------------
    y_raw = y_raw - mean(y_raw);
    sy = std(y_raw);
    
    if sy == 0
        error('Standard deviation of y_raw is zero.');
    end
    
    y = y_raw / sy;
    beta_true_std = beta_true / sy;
    
    % -----------------------------
    % Store outputs
    % -----------------------------
    summary_list = {};
    
    % ============================================================
    % Method runner (internal function handle)
    % ============================================================
    run_method = @(prior_name) run_bayesreg_method( ...
        X, y, beta_true_std, burn, iter, thin, ...
        use_topk, tol, prior_name);
    
    % ============================================================
    % Run Bayesian Lasso
    % ============================================================
    if run_bayeslasso == 1
        disp('Running Bayesian Lasso...');
        out = run_method('lasso');
        out.Method = "Bayesian_Lasso";
        summary_list{end+1} = out;
    end
    
    % ============================================================
    % Run Horseshoe
    % ============================================================
    if run_horseshoe == 1
        disp('Running Horseshoe...');
        out = run_method('horseshoe');
        out.Method = "Horseshoe";
        summary_list{end+1} = out;
    end
    
    % ============================================================
    % Run Horseshoe+
    % ============================================================
    if run_horseshoeplus == 1
        disp('Running Horseshoe+...');
        out = run_method('horseshoe+');
        out.Method = "HorseshoePlus";
        summary_list{end+1} = out;
    end
    
    % ============================================================
    % Add metadata and save
    % ============================================================
    for i = 1:length(summary_list)
        summary_list{i}.Scenario = scenario;
        summary_list{i}.rep = rep;
        summary_list{i}.n = n;
        summary_list{i}.p = p;
    end
    
    % Save combined summary
    summary_tab = struct2table(summary_list{1});
    for i = 2:numel(summary_list)
        summary_tab = [summary_tab; struct2table(summary_list{i})]; %#ok<AGROW>
    end
    
    outfile_summary = sprintf( ...
        '%s/Output_scen_%d_n_%d_p_%d_BayesRegALL_rep_%d.csv', ...
        results_dir, scenario, n, p, rep);
    
    writetable(summary_tab, outfile_summary);
    disp(summary_tab);
    
end

% ============================================================
% ----------- Helper Functions -------------------------------
% ============================================================

function out = run_bayesreg_method(X, y, beta_true, burn, iter, thin, use_topk, tol, prior)

tStart = tic;

[beta_post, ~, retval] = bayesreg( ...
    X, y, ...
    'gaussian', ...
    prior, ...
    'burnin', burn, ...
    'nsamples', iter - burn, ...
    'thin', thin, ...
    'display', false);

comp_time = toc(tStart);

% Posterior mean
beta_est = retval.muB(:);

out = evaluate_fit(beta_est, beta_true, X, comp_time, use_topk, tol);
end


function out = evaluate_fit(beta_est, beta_ref, X, comp_time, use_topk, tol)

beta_est = beta_est(:);
beta_ref = beta_ref(:);

% Estimation metrics
l2_beta   = sqrt(sum((beta_est - beta_ref).^2));
l1_beta   = sum(abs(beta_est - beta_ref));
rmse_beta = sqrt(mean((beta_est - beta_ref).^2));

% Prediction
pred_mse = mean((X * beta_ref - X * beta_est).^2);

% Support recovery
true_supp = (beta_ref ~= 0);
p = length(beta_ref);

if use_topk
    k = sum(true_supp);
    est_supp = false(p,1);
    [~, ord] = sort(abs(beta_est),'descend');
    est_supp(ord(1:k)) = true;
else
    est_supp = abs(beta_est) > tol;
end

TP = sum(true_supp & est_supp);
FP = sum(~true_supp & est_supp);
TN = sum(~true_supp & ~est_supp);
FN = sum(true_supp & ~est_supp);

TPR = TP / max(TP + FN,1);
FPR = FP / max(FP + TN,1);
FDR = FP / max(TP + FP,1);

mcc_num = double(TP)*double(TN) - double(FP)*double(FN);
mcc_den = sqrt(double(TP+FP)*double(TP+FN)*double(TN+FP)*double(TN+FN));

if isnan(mcc_den) || mcc_den <= 0
    MCC = NaN;
else
    MCC = mcc_num / mcc_den;
end

out = struct();
out.L2 = l2_beta;
out.L1 = l1_beta;
out.RMSE = rmse_beta;
out.PredMSE = pred_mse;
out.TPR = TPR;
out.FPR = FPR;
out.MCC = MCC;
out.FDR = FDR;
out.CompTime = comp_time;
out.TP = TP;
out.FP = FP;
out.TN = TN;
out.FN = FN;
out.tol = tol;
end

