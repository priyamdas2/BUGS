function [X, y, beta] = generate_scenario_1_ultrahigh(n, true_p, p, sigma, use_single, block_size)
%GENERATE_SCENARIO_1 Generate one sparse regression dataset.
%
% INPUTS
%   n          : sample size
%   true_p     : number of nonzero beta entries (must be <= 10)
%   p          : total number of predictors
%   sigma      : noise standard deviation
%   use_single : logical; if true, store outputs in single precision
%                (default: false)
%   block_size : number of columns generated at a time
%                (default: min(p, 5000))
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
%   [X, y, beta] = generate_scenario_1(100, 5, 100000, 1, true, 5000);

    if nargin < 5 || isempty(use_single)
        use_single = false;
    end
    if nargin < 6 || isempty(block_size)
        block_size = min(p, 5000);
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
    if block_size < 1
        error('block_size must be positive.');
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

    % Precision control
    if use_single
        out_class = 'single';
        beta_template = single(beta_template);
        sigma_cast = single(sigma);
    else
        out_class = 'double';
        sigma_cast = sigma;
    end

    % Construct beta
    beta = zeros(p, 1, out_class);
    beta(1:true_p) = beta_template(1:true_p);

    % Preallocate design matrix
    X = zeros(n, p, out_class);

    % Initialize signal accumulator so we do not need a full X*beta later
    y_signal = zeros(n, 1, out_class);

    % Generate design matrix blockwise
    for j1 = 1:block_size:p
        j2 = min(j1 + block_size - 1, p);
        idx = j1:j2;
        nb = numel(idx);

        % Generate block
        Xblk = randn(n, nb, out_class);

        % Standardize columns:
        Xblk = Xblk - mean(Xblk, 1);
        col_sd = sqrt(sum(Xblk.^2, 1) / n);

        % Very rare safeguard for zero-variance columns
        bad = (col_sd == 0);
        if any(bad)
            Xtmp = randn(n, sum(bad), out_class);
            Xtmp = Xtmp - mean(Xtmp, 1);
            tmp_sd = sqrt(sum(Xtmp.^2, 1) / n);

            % If a pathological zero-variance case still occurs, fix by setting sd to 1
            tmp_sd(tmp_sd == 0) = 1;

            Xblk(:, bad) = Xtmp ./ tmp_sd;
            col_sd(bad) = 1;
        end

        Xblk = Xblk ./ col_sd;

        % Store block
        X(:, idx) = Xblk;

        % Update signal only for active coordinates in this block
        if j1 <= true_p
            active_end = min(j2, true_p);
            active_local = 1:(active_end - j1 + 1);
            y_signal = y_signal + Xblk(:, active_local) * beta(idx(active_local));
        end
    end

    % Generate response
    eps = sigma_cast * randn(n, 1, out_class);
    y = y_signal + eps;
end