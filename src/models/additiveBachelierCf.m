function phi = additiveBachelierCf(u, T, sigmaT, eta, k)
% The function implements the cf of the BA model for alfa=1/2

% Important: sigma_t is not the ATM sigma

psi = @(u,kappa,alfa) (1./kappa).*((1-alfa)./alfa) .* (1 - (1 + (u.*kappa)./(1-alfa)).^alfa);

log_psi = psi(1i.*u.*eta.*sigmaT.*sqrt(T) + 0.5.*(u.^2).*(sigmaT.^2).*T,k,0.5) + 1i.*u.*eta.*sigmaT.*sqrt(T);

phi = exp(log_psi);
end