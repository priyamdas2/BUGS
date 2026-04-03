function metrics = compute_metrics_from_beta_BUGS(X, beta_true, beta_hat, comp_time, sel_thresh)
% Compute estimation, prediction, selection, and timing metrics.

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

    % Support recovery
    true_support = (beta_true ~= 0);
    est_support  = (abs(beta_hat) > sel_thresh);

    TP = sum(est_support & true_support);
    FP = sum(est_support & ~true_support);
    TN = sum(~est_support & ~true_support);
    FN = sum(~est_support & true_support);

    if (TP + FN) > 0
        TPR = TP / (TP + FN);
    else
        TPR = NaN;
    end

    if (FP + TN) > 0
        FPR = FP / (FP + TN);
    else
        FPR = NaN;
    end

    if (TP + FP) > 0
        FDR = FP / (TP + FP);
    else
        FDR = 0;
    end

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

    % Optional counts
    metrics.TP = TP;
    metrics.FP = FP;
    metrics.TN = TN;
    metrics.FN = FN;
    metrics.sel_thresh = sel_thresh;
end