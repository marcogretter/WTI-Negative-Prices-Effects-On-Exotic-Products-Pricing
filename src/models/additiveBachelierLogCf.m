function logphi = additiveBachelierLogCf(u, T, sigmaT, eta, k)
% additiveBachelierLogCf
%
% Log-characteristic function of the additive Bachelier model
% for alpha = 1/2.
%
% Same convention as additiveBachelierCf, but returns log(phi).

    alpha = 0.5;

    psi = @(z,kappa,alfa) (1./kappa) .* ((1-alfa)./alfa) .* (1 - (1 + ((z.*kappa)./(1-alfa))).^alfa);

    z = 1i .* u .* eta .* sigmaT .* sqrt(T) + 0.5 .* (u.^2) .* (sigmaT.^2) .* T;

    logphi = psi(z, k, alpha) + 1i .* u .* eta .* sigmaT .* sqrt(T);

end