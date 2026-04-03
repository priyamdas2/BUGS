function [X, y, beta] = generate_scenario_2(n, true_p, p, sigma, rho)
%GENERATE_SCENARIO_2 Generate one sparse regression dataset with
%Toeplitz-correlated predictors.
%
% INPUTS
%   n       : sample size
%   true_p  : number of nonzero beta entries (must be <= 10)
%   p       : total number of predictors
%   sigma   : noise standard deviation
%   rho     : Toeplitz correlation parameter, must satisfy |rho| < 1
%
% OUTPUTS
%   X     : n x p design matrix (standardized columns)
%   y     : n x 1 response vector
%   beta  : p x 1 true coefficient vector
%
% DETAILS
%   - Rows of X are generated from a Gaussian Toeplitz covariance
%         Sigma(j,k) = rho^|j-k|
%     using a fast AR(1)-style recursion across columns, so this scales
%     well to very large p without forming the full p x p covariance.
%
%   - After generation, each column of X is standardized so that:
%         sum_i x_ij = 0
%         (1/n) sum_i x_ij^2 = 1
%
%   - beta has first true_p entries nonzero, using a moderate decaying
%     template intended to create a meaningful correlated sparse-regression
%     benchmark.
%
%   - y = X * beta + eps, eps ~ N(0, sigma^2 I)
%
% EXAMPLE
%   [X, y, beta] = generate_scenario_2(100, 5, 200, 1, 0.5);

    if true_p > 10
        error('true_p must be <= 10.');
    end
    if true_p > p
        error('true_p must be <= p.');
    end
    if n < 2
        error('n must be at least 2.');
    end
    if sigma <= 0
        error('sigma must be positive.');
    end
    if abs(rho) >= 1
        error('rho must satisfy abs(rho) < 1.');
    end

    % Moderate decaying signal template for correlated design
    beta_template = [ ...
         2.5; ...
        -2.0; ...
         1.8; ...
        -1.5; ...
         1.2; ...
         1.0; ...
        -0.9; ...
         0.8; ...
         0.7; ...
        -0.7 ];

    % ---------------------------------------------------------
    % Generate Toeplitz-correlated design matrix
    % Each row follows an AR(1)-style recursion across columns:
    %   X(i,1) ~ N(0,1)
    %   X(i,j) = rho * X(i,j-1) + sqrt(1-rho^2) * eps_ij
    % This yields Cov(X_ij, X_ik) = rho^|j-k|.
    % ---------------------------------------------------------
    X = zeros(n, p);
    X(:,1) = randn(n, 1);

    if p > 1
        innovation_sd = sqrt(1 - rho^2);
        for j = 2:p
            X(:,j) = rho * X(:,j-1) + innovation_sd * randn(n, 1);
        end
    end

    % ---------------------------------------------------------
    % Standardize columns:
    %   mean zero and empirical variance one, i.e.
    %   (1/n) sum_i x_ij^2 = 1
    % ---------------------------------------------------------
    X = X - mean(X, 1);
    col_sd = sqrt(sum(X.^2, 1) / n);

    % Safeguard for zero-variance columns
    bad = (col_sd == 0);
    while any(bad)
        nb = sum(bad);

        % Regenerate bad columns independently, then standardize again.
        X(:, bad) = randn(n, nb);
        X(:, bad) = X(:, bad) - mean(X(:, bad), 1);
        col_sd(bad) = sqrt(sum(X(:, bad).^2, 1) / n);

        bad = (col_sd == 0);
    end

    X = X ./ col_sd;

    % ---------------------------------------------------------
    % Construct beta
    % ---------------------------------------------------------
    beta = zeros(p, 1);
    beta(1:true_p) = beta_template(1:true_p);

    % ---------------------------------------------------------
    % Generate response
    % ---------------------------------------------------------
    eps = sigma * randn(n, 1);
    y = X * beta + eps;
end