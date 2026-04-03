function out = BUGS_Active_mcmc(X, y, opts)
% BUGS-Active: sparsity-aware approximate MCMC for
% Bayesian Univariate-Guided Shrinkage (BUGS / B-UGS)
%
% IMPORTANT:
% This is a Johndrow-inspired active-set approximation, not an exact
% implementation of Johndrow et al. (2020).
%
% Main idea:
%   - beta update remains global
%   - lambda_j updates are restricted to an active set
%   - inactive lambda_j are reset to a small baseline value
%   - full beta draws can be disabled for memory efficiency
%
% INPUTS
%   X      : n x p standardized design
%   y      : n x 1 standardized response
%   opts   : struct with fields
%       % shared with BUGS_mcmc
%       .niter
%       .burnin
%       .thin
%       .tau0
%       .a_c, .b_c
%       .a_sigma, .b_sigma
%       .sigma_eta
%       .eps_score
%       .Cz
%       .slice_w
%       .slice_m
%       .seed
%       .verbose
%       .save_lambda_max
%
%       % active-set specific
%       .active_beta_thresh
%       .active_cap
%       .min_active
%       .guidance_keep
%       .lambda2_inactive
%
%       % memory/saving specific
%       .save_full_beta
%       .max_p_full_save
%       .save_top_beta
%
% OUTPUT
%   out struct with posterior draws and diagnostics

if nargin < 3
    opts = struct();
end
opts = set_defaults_active(opts);

rng(opts.seed);

[n, p] = size(X);
y = y(:);

%--------------------------------------------------------------
% Compute guidance statistics zstar
%--------------------------------------------------------------
s = abs(X' * y) / n;
ztilde = log(s + opts.eps_score);
z = (ztilde - mean(ztilde)) / std(ztilde);
zstar = max(-opts.Cz, min(z, opts.Cz));
zstar = zstar(:);

%--------------------------------------------------------------
% Initialize parameters
%--------------------------------------------------------------
beta = zeros(p,1);
lambda2 = ones(p,1);
tau2 = opts.tau0^2 / 10;
c2 = 1.0;
eta = 0.1;
sigma2 = 1.0;

In = speye(n);

%--------------------------------------------------------------
% Saving setup
%--------------------------------------------------------------
nkeep = floor((opts.niter - opts.burnin) / opts.thin);

save_full_beta = opts.save_full_beta && (p <= opts.max_p_full_save);
if save_full_beta
    beta_save = zeros(nkeep, p);
else
    beta_save = [];
end

if opts.save_top_beta > 0
    beta_top_save = zeros(nkeep, opts.save_top_beta);
    beta_top_idx_save = zeros(nkeep, opts.save_top_beta);
else
    beta_top_save = [];
    beta_top_idx_save = [];
end

sigma2_save = zeros(nkeep, 1);
tau_save    = zeros(nkeep, 1);
c2_save     = zeros(nkeep, 1);
eta_save    = zeros(nkeep, 1);
active_size_save = zeros(nkeep, 1);
lambda_save = zeros(nkeep, min(p, opts.save_lambda_max));

keep_idx = 0;

if opts.verbose
    fprintf('Running %d iterations (BUGS-Active)...\n', opts.niter);
end

for iter = 1:opts.niter

    %----------------------------------------------------------
    % Update kappa_tilde_j^2
    %----------------------------------------------------------
    exp_term = exp(eta * zstar);
    numer = c2 * tau2 .* lambda2 .* exp_term;
    kappa2 = numer ./ (c2 + tau2 .* lambda2 .* exp_term);
    kappa2 = max(kappa2, 1e-12);

    %----------------------------------------------------------
    % Update beta
    %----------------------------------------------------------
    beta = sample_beta_woodbury(X, y, sigma2, kappa2, In);

    %----------------------------------------------------------
    % Update sigma2
    %----------------------------------------------------------
    resid = y - X * beta;
    shape_sig = opts.a_sigma + 0.5 * (n + p);
    rate_sig  = opts.b_sigma + 0.5 * (resid' * resid) + 0.5 * sum(beta.^2 ./ kappa2);
    sigma2 = 1 / gamrnd(shape_sig, 1 / rate_sig);

    %----------------------------------------------------------
    % Construct active set
    %
    % A coordinate is active if it currently looks important enough
    % to receive a full lambda_j update.
    %----------------------------------------------------------
    absbeta = abs(beta);

    % 1) threshold by current |beta|
    active_mask = absbeta > opts.active_beta_thresh;

    % 2) restrict to top active_cap by |beta| if requested
    if opts.active_cap < p
        [~, ord_beta] = sort(absbeta, 'descend');
        top_idx = ord_beta(1:min(opts.active_cap, p));
        cap_mask = false(p,1);
        cap_mask(top_idx) = true;
        active_mask = active_mask & cap_mask;
    end

    % 3) always include top guidance coordinates if requested
    if opts.guidance_keep > 0
        [~, ord_z] = sort(abs(zstar), 'descend');
        keep_z_idx = ord_z(1:min(opts.guidance_keep, p));
        active_mask(keep_z_idx) = true;
    end

    % 4) enforce minimum active size
    if opts.min_active > 0 && nnz(active_mask) < opts.min_active
        [~, ord_beta] = sort(absbeta, 'descend');
        active_mask(ord_beta(1:min(opts.min_active, p))) = true;
    end

    active_idx = find(active_mask);
    inactive_idx = find(~active_mask);

    %----------------------------------------------------------
    % Update local lambda_j^2 only on active set
    %----------------------------------------------------------
    for kk = 1:numel(active_idx)
        j = active_idx(kk);
        theta0 = log(lambda2(j));
        logf = @(th) logpost_theta_j(th, beta(j), sigma2, tau2, c2, eta, zstar(j));
        theta_new = slice_sample_univariate(theta0, logf, opts.slice_w, opts.slice_m, -Inf, Inf);
        lambda2(j) = exp(theta_new);
    end

    %----------------------------------------------------------
    % Reset inactive lambda_j^2 to small baseline
    %----------------------------------------------------------
    lambda2(inactive_idx) = opts.lambda2_inactive;

    %----------------------------------------------------------
    % Update tau
    %----------------------------------------------------------
    phi0 = log(tau2);
    logf_phi = @(phi) logpost_phi(phi, beta, sigma2, lambda2, c2, eta, zstar, opts.tau0);
    phi_new = slice_sample_univariate(phi0, logf_phi, opts.slice_w, opts.slice_m, -Inf, Inf);
    tau2 = exp(phi_new);

    %----------------------------------------------------------
    % Update c2
    %----------------------------------------------------------
    psi0 = log(c2);
    logf_psi = @(psi) logpost_psi(psi, beta, sigma2, lambda2, tau2, eta, zstar, opts.a_c, opts.b_c);
    psi_new = slice_sample_univariate(psi0, logf_psi, opts.slice_w, opts.slice_m, -Inf, Inf);
    c2 = exp(psi_new);

    %----------------------------------------------------------
    % Update eta
    %----------------------------------------------------------
    logf_eta = @(et) logpost_eta(et, beta, sigma2, lambda2, tau2, c2, zstar, opts.sigma_eta);
    eta = slice_sample_univariate(eta, logf_eta, opts.slice_w, opts.slice_m, 0, Inf);

    %----------------------------------------------------------
    % Save
    %----------------------------------------------------------
    if iter > opts.burnin && mod(iter - opts.burnin, opts.thin) == 0
        keep_idx = keep_idx + 1;

        if save_full_beta
            beta_save(keep_idx,:) = beta';
        end

        if opts.save_top_beta > 0
            [~, ord] = sort(abs(beta), 'descend');
            top_idx = ord(1:min(opts.save_top_beta, p));
            beta_top_idx_save(keep_idx,1:numel(top_idx)) = top_idx';
            beta_top_save(keep_idx,1:numel(top_idx)) = beta(top_idx)';
        end

        sigma2_save(keep_idx) = sigma2;
        tau_save(keep_idx) = sqrt(tau2);
        c2_save(keep_idx) = c2;
        eta_save(keep_idx) = eta;
        active_size_save(keep_idx) = numel(active_idx);
        lambda_save(keep_idx,:) = sqrt(lambda2(1:min(p, opts.save_lambda_max)))';
    end

    if opts.verbose && (mod(iter, opts.display_every) == 0)
        fprintf('Iter %d/%d | sigma2=%.4f tau=%.4f c2=%.4f eta=%.4f | active=%d\n', ...
            iter, opts.niter, sigma2, sqrt(tau2), c2, eta, numel(active_idx));
    end
end

%--------------------------------------------------------------
% Posterior summaries
%--------------------------------------------------------------
out.beta_draws = beta_save;

if save_full_beta
    out.beta_mean = mean(beta_save, 1)';
    out.beta_median = median(beta_save, 1)';
else
    out.beta_mean = [];
    out.beta_median = [];
end

out.beta_top_draws = beta_top_save;
out.beta_top_idx_draws = beta_top_idx_save;

out.sigma2_draws = sigma2_save;
out.tau_draws = tau_save;
out.c2_draws = c2_save;
out.eta_draws = eta_save;
out.active_size_draws = active_size_save;
out.lambda_draws_partial = lambda_save;
out.zstar = zstar;
out.opts = opts;

end

%======================================================================
% Helpers
%======================================================================

function opts = set_defaults_active(opts)

%--------------------------------------------------------------
% Defaults shared with BUGS_mcmc: EXACTLY MATCH
%--------------------------------------------------------------
def.niter = 6000;
def.burnin = 3000;
def.thin = 5;
def.tau0 = 0.02;

def.a_c = 2;
def.b_c = 2;
def.a_sigma = 2;
def.b_sigma = 2;
def.sigma_eta = 0.1;

def.eps_score = 1e-6;
def.Cz = 3;
def.slice_w = 0.5;
def.slice_m = 100;

def.seed = 1;
def.verbose = true;
def.save_lambda_max = 100;
def.display_every = 100;

%--------------------------------------------------------------
% Active-set specific defaults
%--------------------------------------------------------------
def.active_beta_thresh = 0.01;
def.active_cap = 500;
def.min_active = 25;
def.guidance_keep = 50;
def.lambda2_inactive = 1e-6;

%--------------------------------------------------------------
% Memory / saving defaults
%--------------------------------------------------------------
def.save_full_beta = false;
def.max_p_full_save = 5000;
def.save_top_beta = 100;

fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i})
        opts.(fn{i}) = def.(fn{i});
    end
end
end

function beta = sample_beta_woodbury(X, y, sigma2, kappa2, In)
% Exact Gaussian sampler for
% beta | - ~ N(m, sigma2 * (X'X + D^{-1})^{-1}),
% where D = diag(kappa2)

D = kappa2(:);
sigma = sqrt(sigma2);

u = randn(length(D),1) .* sqrt(D);
delta = randn(size(X,1),1);

XD = X .* (D');
A = In + XD * X';
rhs = y / sigma - X * u - delta;
w = A \ rhs;

beta = sigma * (u + D .* (X' * w));
end

function val = logpost_theta_j(theta, beta_j, sigma2, tau2, c2, eta, zstar_j)
lambda2 = exp(theta);
exp_term = exp(eta * zstar_j);
kappa2 = (c2 * tau2 * lambda2 * exp_term) / (c2 + tau2 * lambda2 * exp_term);
kappa2 = max(kappa2, 1e-12);

val = -0.5 * log(kappa2) ...
    - (beta_j^2) / (2 * sigma2 * kappa2) ...
    + 0.5 * theta - log1pexp(theta);
end

function val = logpost_phi(phi, beta, sigma2, lambda2, c2, eta, zstar, tau0)
tau2 = exp(phi);
exp_term = exp(eta * zstar(:));
kappa2 = (c2 * tau2 .* lambda2(:) .* exp_term) ./ ...
    (c2 + tau2 .* lambda2(:) .* exp_term);
kappa2 = max(kappa2, 1e-12);

logprior = -log(1 + tau2 / tau0^2) + 0.5 * phi;

val = -0.5 * sum(log(kappa2)) ...
    - 0.5 / sigma2 * sum(beta(:).^2 ./ kappa2) ...
    + logprior;
end

function val = logpost_psi(psi, beta, sigma2, lambda2, tau2, eta, zstar, a_c, b_c)
c2 = exp(psi);
exp_term = exp(eta * zstar(:));
kappa2 = (c2 * tau2 .* lambda2(:) .* exp_term) ./ ...
    (c2 + tau2 .* lambda2(:) .* exp_term);
kappa2 = max(kappa2, 1e-12);

logprior = -(a_c + 1) * psi - b_c / c2 + psi;

val = -0.5 * sum(log(kappa2)) ...
    - 0.5 / sigma2 * sum(beta(:).^2 ./ kappa2) ...
    + logprior;
end

function val = logpost_eta(eta, beta, sigma2, lambda2, tau2, c2, zstar, sigma_eta)
if eta < 0
    val = -Inf;
    return;
end

exp_term = exp(eta * zstar(:));
kappa2 = (c2 * tau2 .* lambda2(:) .* exp_term) ./ ...
    (c2 + tau2 .* lambda2(:) .* exp_term);
kappa2 = max(kappa2, 1e-12);

logprior = -(eta^2) / (2 * sigma_eta^2);

val = -0.5 * sum(log(kappa2)) ...
    - 0.5 / sigma2 * sum(beta(:).^2 ./ kappa2) ...
    + logprior;
end

function x_new = slice_sample_univariate(x0, logf, w, m, Lbound, Ubound)
logy = logf(x0) - exprnd(1);

u = rand();
L = x0 - w * u;
R = L + w;

if isfinite(Lbound), L = max(L, Lbound); end
if isfinite(Ubound), R = min(R, Ubound); end

J = floor(m * rand());
K = (m - 1) - J;

while J > 0 && L > Lbound && logf(L) > logy
    L = L - w;
    if isfinite(Lbound), L = max(L, Lbound); end
    J = J - 1;
end
while K > 0 && R < Ubound && logf(R) > logy
    R = R + w;
    if isfinite(Ubound), R = min(R, Ubound); end
    K = K - 1;
end

while true
    x1 = L + rand() * (R - L);
    lf = logf(x1);
    if isfinite(lf) && lf >= logy
        x_new = x1;
        return;
    else
        if x1 < x0
            L = x1;
        else
            R = x1;
        end
    end
end
end

function y = log1pexp(x)
y = zeros(size(x));
idx = x > 0;
y(idx) = x(idx) + log1p(exp(-x(idx)));
y(~idx) = log1p(exp(x(~idx)));
end