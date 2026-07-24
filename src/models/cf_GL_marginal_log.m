function logphi = cf_GL_marginal_log(u, t, sigma_t, alpha, beta)
% cf_GL_marginal_log
%
% Log marginal CF of f_t = sigma_t sqrt(t) Z_GL.

    q = sigma_t * sqrt(t);

    logphi = generalizedLogisticLogCf(q .* u, alpha, beta);

end