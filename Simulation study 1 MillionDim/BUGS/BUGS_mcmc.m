function out = BUGS_mcmc(X, y, opts)
% Bayesian Univariate-Guided Shrinkage (BUGS / B-UGS) MCMC
%
% Model:
%   y = X*beta + eps,  eps ~ N(0, sigma2 I)
%   beta_j | ... ~ N(0, sigma2 * kappa_tilde_j^2)
%   kappa_tilde_j^2 = [c2 * tau2 * lambda2_j * exp(eta*zstar_j)] ...
%                     / [c2 + tau2 * lambda2_j * exp(eta*zstar_j)]
%
% Priors:
%   lambda_j ~ half-Cauchy(0,1)
%   tau      ~ half-Cauchy(0,tau0)
%   c2       ~ IG(a_c, b_c)
%   eta      ~ half-Normal(0, sigma_eta^2)
%   sigma2   ~ IG(a_sigma, b_sigma)
%
% INPUTS
%   X      : n x p standardized design
%   y      : n x 1 standardized response
%   opts   : struct with fields
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
%
% OUTPUT
%   out struct with posterior draws and diagnostics

if nargin < 3
    opts = struct();
end
opts = set_defaults(opts);

rng(opts.seed);

[n, p] = size(X);
y = y(:);

%--------------------------------------------------------------
% Compute guidance statistics zstar from the draft
% s_j = |x_j' y| / n
% ztilde_j = log(s_j + eps_score)
% z_j standardized and clipped to [-Cz, Cz]
%--------------------------------------------------------------
s = abs(X' * y) / n;
ztilde = log(s + opts.eps_score);
z = (ztilde - mean(ztilde)) / std(ztilde);
zstar = max(-opts.Cz, min(z, opts.Cz));

%--------------------------------------------------------------
% Initialize parameters
%--------------------------------------------------------------
beta = zeros(p,1);
lambda2 = ones(p,1);          % lambda_j^2
tau2 = opts.tau0^2 / 10;      % initialize modestly
c2 = 1.0;
eta = 0.1;
sigma2 = 1.0;

% Precompute
Xt = X';
In = speye(n);

nkeep = floor((opts.niter - opts.burnin) / opts.thin);
beta_save   = zeros(nkeep, p);
sigma2_save = zeros(nkeep, 1);
tau_save    = zeros(nkeep, 1);
c2_save     = zeros(nkeep, 1);
eta_save    = zeros(nkeep, 1);
lambda_save = zeros(nkeep, min(p, opts.save_lambda_max)); % partial save if p huge

keep_idx = 0;

if opts.verbose
    fprintf('Running %d iterations...\n', opts.niter);
end

for iter = 1:opts.niter
    
    %----------------------------------------------------------
    % Update kappa_tilde_j^2 and D
    %----------------------------------------------------------
    exp_term = exp(eta * zstar(:));
    numer = c2 * tau2 .* lambda2 .* exp_term;
    kappa2 = numer ./ (c2 + tau2 .* lambda2 .* exp_term);
    kappa2 = max(kappa2, 1e-12); % numerical safety
    
    %----------------------------------------------------------
    % Update beta using Woodbury-based sampler
    % beta | - ~ N(m_beta, sigma2 * (X'X + D^{-1})^{-1})
    %----------------------------------------------------------
    beta = sample_beta_woodbury(X, y, sigma2, kappa2, In);
    
    %----------------------------------------------------------
    % Update sigma2 from IG
    % shape = a_sigma + (n+p)/2
    % rate  = b_sigma + 0.5||y-Xb||^2 + 0.5 sum beta_j^2 / kappa2_j
    %----------------------------------------------------------
    resid = y - X * beta;
    shape_sig = opts.a_sigma + 0.5 * (n + p);
    rate_sig  = opts.b_sigma + 0.5 * (resid' * resid) + 0.5 * sum(beta.^2 ./ kappa2);
    sigma2 = 1 / gamrnd(shape_sig, 1 / rate_sig);
    
    %----------------------------------------------------------
    % Update local lambda_j^2 via slice on theta_j = log(lambda_j^2)
    % log p(theta_j | -) =
    %   -0.5 log(kappa2_j) - beta_j^2/(2 sigma2 kappa2_j)
    %   + 0.5 theta_j - log(1 + exp(theta_j))
    %----------------------------------------------------------
    for j = 1:p
        theta0 = log(lambda2(j));
        logf = @(th) logpost_theta_j(th, beta(j), sigma2, tau2, c2, eta, zstar(j));
        theta_new = slice_sample_univariate(theta0, logf, opts.slice_w, opts.slice_m, -Inf, Inf);
        lambda2(j) = exp(theta_new);
    end
    
    %----------------------------------------------------------
    % Update tau via slice on phi = log(tau^2)
    % Include half-Cauchy(0, tau0) prior for tau
    % with Jacobian from phi -> tau
    %----------------------------------------------------------
    phi0 = log(tau2);
    logf_phi = @(phi) logpost_phi(phi, beta, sigma2, lambda2, c2, eta, zstar, opts.tau0);
    phi_new = slice_sample_univariate(phi0, logf_phi, opts.slice_w, opts.slice_m, -Inf, Inf);
    tau2 = exp(phi_new);
    
    %----------------------------------------------------------
    % Update c2 via slice on psi = log(c2)
    % Include IG(a_c, b_c) prior and Jacobian
    %----------------------------------------------------------
    psi0 = log(c2);
    logf_psi = @(psi) logpost_psi(psi, beta, sigma2, lambda2, tau2, eta, zstar, opts.a_c, opts.b_c);
    psi_new = slice_sample_univariate(psi0, logf_psi, opts.slice_w, opts.slice_m, -Inf, Inf);
    c2 = exp(psi_new);
    
    %----------------------------------------------------------
    % Update eta via slice on [0, Inf)
    % half-Normal prior: log pi(eta) = const - eta^2/(2 sigma_eta^2), eta >= 0
    %----------------------------------------------------------
    logf_eta = @(et) logpost_eta(et, beta, sigma2, lambda2, tau2, c2, zstar, opts.sigma_eta);
    eta = slice_sample_univariate(eta, logf_eta, opts.slice_w, opts.slice_m, 0, Inf);
    
    %----------------------------------------------------------
    % Save
    %----------------------------------------------------------
    if iter > opts.burnin && mod(iter - opts.burnin, opts.thin) == 0
        keep_idx = keep_idx + 1;
        beta_save(keep_idx,:) = beta';
        sigma2_save(keep_idx) = sigma2;
        tau_save(keep_idx) = sqrt(tau2);
        c2_save(keep_idx) = c2;
        eta_save(keep_idx) = eta;
        lambda_save(keep_idx,:) = sqrt(lambda2(1:min(p, opts.save_lambda_max)))';
    end
    
    if opts.verbose && (mod(iter, opts.display_every) == 0)
        fprintf('Iter %d/%d | sigma2=%.4f tau=%.4f c2=%.4f eta=%.4f\n', ...
            iter, opts.niter, sigma2, sqrt(tau2), c2, eta);
    end
end

% Posterior summaries
out.beta_draws = beta_save;
out.beta_mean = mean(beta_save, 1)';
out.beta_median = median(beta_save, 1)';
out.sigma2_draws = sigma2_save;
out.tau_draws = tau_save;
out.c2_draws = c2_save;
out.eta_draws = eta_save;
out.lambda_draws_partial = lambda_save;
out.zstar = zstar;
out.opts = opts;
end

%======================================================================
% Helpers
%======================================================================

function opts = set_defaults(opts)
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

    % Draw auxiliary variables
    u = randn(length(D),1) .* sqrt(D);   % u ~ N(0, D)
    delta = randn(size(X,1),1);          % delta ~ N(0, I_n)

    % Solve (I + X D X') w = y/sigma - X u - delta
    XD = X .* (D');
    A = In + XD * X';
    rhs = y / sigma - X * u - delta;
    w = A \ rhs;

    % Final draw
    beta = sigma * (u + D .* (X' * w));
end

function val = logpost_theta_j(theta, beta_j, sigma2, tau2, c2, eta, zstar_j)
% theta = log(lambda_j^2)
lambda2 = exp(theta);
exp_term = exp(eta * zstar_j);
kappa2 = (c2 * tau2 * lambda2 * exp_term) / (c2 + tau2 * lambda2 * exp_term);
kappa2 = max(kappa2, 1e-12);

% half-Cauchy on lambda -> density for lambda2 up to const:
% p(theta) proportional exp(theta/2)/(1+exp(theta))
val = -0.5 * log(kappa2) ...
    - (beta_j^2) / (2 * sigma2 * kappa2) ...
    + 0.5 * theta - log1pexp(theta);
end

function val = logpost_phi(phi, beta, sigma2, lambda2, c2, eta, zstar, tau0)
% phi = log(tau^2), tau = exp(phi/2)
tau2 = exp(phi);
exp_term = exp(eta * zstar(:));
kappa2 = (c2 * tau2 .* lambda2(:) .* exp_term) ./ ...
    (c2 + tau2 .* lambda2(:) .* exp_term);
kappa2 = max(kappa2, 1e-12);

% half-Cauchy prior on tau with scale tau0
tau = sqrt(tau2);
% log p(tau) + log|d tau / d phi| = -log(1 + (tau/tau0)^2) + 0.5 phi + const
logprior = -log(1 + tau2 / tau0^2) + 0.5 * phi;

val = -0.5 * sum(log(kappa2)) ...
    - 0.5 / sigma2 * sum(beta(:).^2 ./ kappa2) ...
    + logprior;
end

function val = logpost_psi(psi, beta, sigma2, lambda2, tau2, eta, zstar, a_c, b_c)
% psi = log(c2)
c2 = exp(psi);
exp_term = exp(eta * zstar(:));
kappa2 = (c2 * tau2 .* lambda2(:) .* exp_term) ./ ...
    (c2 + tau2 .* lambda2(:) .* exp_term);
kappa2 = max(kappa2, 1e-12);

% IG(a_c, b_c) on c2: p(c2) proportional c2^(-a_c-1) exp(-b_c/c2)
% plus Jacobian +psi
logprior = -(a_c + 1) * psi - b_c / c2 + psi;

val = -0.5 * sum(log(kappa2)) ...
    - 0.5 / sigma2 * sum(beta(:).^2 ./ kappa2) ...
    + logprior;
end

function val = logpost_eta(eta, beta, sigma2, lambda2, tau2, c2, zstar, sigma_eta)
% eta >= 0
if eta < 0
    val = -Inf;
    return;
end

exp_term = exp(eta * zstar(:));
kappa2 = (c2 * tau2 .* lambda2(:) .* exp_term) ./ ...
    (c2 + tau2 .* lambda2(:) .* exp_term);
kappa2 = max(kappa2, 1e-12);

% half-normal prior
logprior = -(eta^2) / (2 * sigma_eta^2);

val = -0.5 * sum(log(kappa2)) ...
    - 0.5 / sigma2 * sum(beta(:).^2 ./ kappa2) ...
    + logprior;
end

function x_new = slice_sample_univariate(x0, logf, w, m, Lbound, Ubound)
% Generic stepping-out slice sampler for 1D target proportional exp(logf(x))
%
% x0     current state
% logf   handle
% w      step size
% m      max stepping-out steps
% bounds lower/upper support bounds (use -Inf/Inf if unbounded)

logy = logf(x0) - exprnd(1);

% Step out
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

% Shrinkage
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
% Numerically stable log(1+exp(x))
y = zeros(size(x));
idx = x > 0;
y(idx) = x(idx) + log1p(exp(-x(idx)));
y(~idx) = log1p(exp(x(~idx)));
end