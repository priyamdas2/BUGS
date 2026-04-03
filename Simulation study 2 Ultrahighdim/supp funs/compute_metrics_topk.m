function metrics = compute_metrics_topk(X, beta_true, beta_hat, comp_time, k)
% Compute L2, L1, RMSE, PredMSE, TPR, FPR, MCC, FDR, CompTime
% using top-k selection based on |beta_hat|.

    beta_true = beta_true(:);
    beta_hat  = beta_hat(:);

    err = beta_hat - beta_true;

    % Estimation metrics
    L2_beta   = norm(err, 2);
    L1_beta   = norm(err, 1);
    RMSE_beta = sqrt(mean(err.^2));

    % Prediction MSE
    pred_err = X * err;
    PredMSE  = mean(pred_err.^2);

    % True support
    true_support = (beta_true ~= 0);

    % Estimated support: top-k by absolute value
    p = length(beta_hat);
    est_support = false(p,1);
    [~, idx] = sort(abs(beta_hat), 'descend');
    est_support(idx(1:k)) = true;

    TP = sum(est_support & true_support);
    FP = sum(est_support & ~true_support);
    TN = sum(~est_support & ~true_support);
    FN = sum(~est_support & true_support);

    TPR = TP / (TP + FN);
    FPR = FP / (FP + TN);
    FDR = FP / (TP + FP);

    denom = sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    if denom > 0
        MCC = (TP*TN - FP*FN) / denom;
    else
        MCC = NaN;
    end

    metrics = struct();
    metrics.L2_beta   = L2_beta;
    metrics.L1_beta   = L1_beta;
    metrics.RMSE_beta = RMSE_beta;
    metrics.PredMSE   = PredMSE;
    metrics.TPR       = TPR;
    metrics.FPR       = FPR;
    metrics.MCC       = MCC;
    metrics.FDR       = FDR;
    metrics.CompTime  = comp_time;
end