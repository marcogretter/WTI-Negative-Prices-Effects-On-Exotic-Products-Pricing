function phi = cf_MA_marginal(u, t, sigma_t, alpha, beta)
% Marginal CF of f_t in Minimal Additive model.
%
% f_t = sigma_t * sqrt(t) * zeta
% zeta has asymmetric Laplace distribution with E[zeta]=0.

    q = sigma_t * sqrt(t);

    gammaLoc = 1/alpha - 1/beta;

    v = q .* u;

    phi = exp(1i .* gammaLoc .* v) ./ (1 - 1i .* v .* (1/beta - 1/alpha) + (v.^2) ./ (alpha .* beta));

end
