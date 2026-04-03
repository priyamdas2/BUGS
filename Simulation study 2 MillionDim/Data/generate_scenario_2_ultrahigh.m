function [X, y, beta] = generate_scenario_2_ultrahigh(n, true_p, p, sigma, rho, use_single)
%GENERATE_SCENARIO_2 Generate one sparse regression dataset with
%Toeplitz-correlated predictors.
%
% INPUTS
%   n          : sample size
%   true_p     : number of nonzero beta entries (must be <= 10)
%   p          : total number of predictors
%   sigma      : noise standard deviation
%   rho        : Toeplitz correlation parameter, must satisfy |rho| < 1
%   use_single : logical; if true, store outputs in single precision
%                (default: false)
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
%   [X, y, beta] = generate_scenario_2(100, 5, 100000, 1, 0.5, true);

    if nargin < 6 || isempty(use_single)
        use_single = false;
    end

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

    % Precision control
    if use_single
        out_class = 'single';
        beta_template = single(beta_template);
        rho_cast = single(rho);
        sigma_cast = single(sigma);
    else
        out_class = 'double';
        rho_cast = rho;
        sigma_cast = sigma;
    end

    % ---------------------------------------------------------
    % Generate Toeplitz-correlated design matrix
    % Each row follows an AR(1)-style recursion across columns:
    %   X(i,1) ~ N(0,1)
    %   X(i,j) = rho * X(i,j-1) + sqrt(1-rho^2) * eps_ij
    % This yields Cov(X_ij, X_ik) = rho^|j-k|.
    %
    % Use filter() row-wise in vectorized form rather than a
    % long explicit MATLAB loop over j = 2,...,p.
    % ---------------------------------------------------------
    innovation_sd = sqrt(1 - rho_cast^2);
    E = randn(n, p, out_class);
    X = filter(1, [1, -rho_cast], E, [], 2);
    X = innovation_sd * X;

    % Replace first column so that X(:,1) ~ N(0,1)
    X(:,1) = randn(n, 1, out_class);

    % ---------------------------------------------------------
    % Standardize columns:
    %   mean zero and empirical variance one, i.e.
    %   (1/n) sum_i x_ij^2 = 1
    % ---------------------------------------------------------
    X = X - mean(X, 1);
    col_sd = sqrt(sum(X.^2, 1) / n);

    % Safeguard for zero-variance columns
    bad = (col_sd == 0);
    if any(bad)
        nb = sum(bad);

        % Regenerate bad columns independently, then standardize again.
        X(:, bad) = randn(n, nb, out_class);
        X(:, bad) = X(:, bad) - mean(X(:, bad), 1);
        col_sd(bad) = sqrt(sum(X(:, bad).^2, 1) / n);

        % Final fallback
        col_sd(col_sd == 0) = 1;
    end

    X = X ./ col_sd;

    % ---------------------------------------------------------
    % Construct beta
    % ---------------------------------------------------------
    beta = zeros(p, 1, out_class);
    beta(1:true_p) = beta_template(1:true_p);

    % ---------------------------------------------------------
    % Generate response
    % ---------------------------------------------------------
    eps = sigma_cast * randn(n, 1, out_class);
    y = X(:,1:true_p) * beta(1:true_p) + eps;
end