function [X, y, beta] = generate_scenario_1(n, true_p, p, sigma)
%GENERATE_SCENARIO_1 Generate one sparse regression dataset.
%
% INPUTS
%   n       : sample size
%   true_p  : number of nonzero beta entries (must be <= 10)
%   p       : total number of predictors
%   sigma   : noise standard deviation
%
% OUTPUTS
%   X     : n x p design matrix (standardized columns)
%   y     : n x 1 response vector
%   beta  : p x 1 true coefficient vector
%
% DETAILS
%   - Each column of X satisfies:
%         sum_i x_ij = 0
%         (1/n) sum_i x_ij^2 = 1
%   - beta has first true_p entries nonzero (fixed template)
%   - y = X * beta + eps, eps ~ N(0, sigma^2 I)
%
% EXAMPLE
%   [X, y, beta] = generate_scenario_1(100, 5, 200, 1);

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

    % Fixed signal template
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

    % Generate design matrix
    X = randn(n, p);

    % Standardize columns:
    X = X - mean(X, 1);
    col_sd = sqrt(sum(X.^2, 1) / n);

    % Safeguard for zero-variance columns
    bad = (col_sd == 0);
    while any(bad)
        X(:, bad) = randn(n, sum(bad));
        X(:, bad) = X(:, bad) - mean(X(:, bad), 1);
        col_sd(bad) = sqrt(sum(X(:, bad).^2, 1) / n);
        bad = (col_sd == 0);
    end

    X = X ./ col_sd;

    % Construct beta
    beta = zeros(p, 1);
    beta(1:true_p) = beta_template(1:true_p);

    % Generate response
    eps = sigma * randn(n, 1);
    y = X * beta + eps;
end