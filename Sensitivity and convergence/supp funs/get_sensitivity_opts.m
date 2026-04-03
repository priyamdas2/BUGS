function [opts, case_label] = get_sensitivity_opts(base_opts, sensitivity_casenum)
%GET_SENSITIVITY_OPTS
% Returns BUGS hyperparameter settings for a given sensitivity-analysis case.
%
% INPUT
%   base_opts            : struct with generic MCMC controls already set
%   sensitivity_casenum  : integer in {1,...,12}
%
% OUTPUT
%   opts       : full options struct for BUGS_mcmc
%   case_label : descriptive text label

opts = base_opts;

% -------------------------------------------------
% Baseline values
% -------------------------------------------------
opts.tau0 = 0.001;
opts.sigma_eta = 0.002;
opts.Cz = 2;
opts.a_c = 10;
opts.b_c = 100;
opts.a_sigma = 2;
opts.b_sigma = 2;

switch sensitivity_casenum

    case 1
        case_label = 'Baseline';

    case 2
        % Stronger global shrinkage
        opts.tau0 = 0.0005;
        case_label = 'Stronger global shrinkage';

    case 3
        % Weaker global shrinkage
        opts.tau0 = 0.005;
        case_label = 'Weaker global shrinkage';

    case 4
        % Guidance nearly off
        opts.sigma_eta = 0.0005;
        case_label = 'Guidance nearly off';

    case 5
        % Moderately freer guidance
        opts.sigma_eta = 0.01;
        case_label = 'Moderately freer guidance';

    case 6
        % Strongly freer guidance
        opts.sigma_eta = 1;
        case_label = 'Strongly freer guidance';

    case 7
        % Strong clipping of z*
        opts.Cz = 1;
        case_label = 'Strong guidance clipping';

    case 8
        % Weak clipping of z*
        opts.Cz = 4;
        case_label = 'Weak guidance clipping';

    case 9
        % Smaller slab mean
        opts.a_c = 10;
        opts.b_c = 50;   % E[c^2] = 50 / 9 ≈ 5.56
        case_label = 'Smaller slab';

    case 10
        % Larger slab mean
        opts.a_c = 10;
        opts.b_c = 500;  % E[c^2] = 500 / 9 ≈ 55.56
        case_label = 'Larger slab';

    case 11
        % Heavier-tailed slab prior
        opts.a_c = 2;
        opts.b_c = 8;    % E[c^2] = 8, but much heavier tail
        case_label = 'Heavier-tailed slab prior';

    case 12
        % Joint permissive / stress-test case
        opts.tau0 = 0.005;
        opts.sigma_eta = 0.05;
        opts.Cz = 4;
        opts.a_c = 10;
        opts.b_c = 500;
        case_label = 'Permissive joint stress test';

    otherwise
        error('Unknown sensitivity_casenum = %d. Must be between 1 and 12.', ...
              sensitivity_casenum);
end
end