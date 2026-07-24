function c = call_MA_normalized(y, alpha, beta)
% Normalized call price for the Minimal Additive model:
%
% c(y; alpha,beta) = E[(Z - y)^+]
%
% where Z has asymmetric Laplace distribution with E[Z] = 0.
%
% alpha = left tail coefficient
% beta  = right tail coefficient

if alpha <= 0 || beta <= 0
    c = NaN(size(y));
    return;
end

y = y(:);

C = 1 / (1/alpha + 1/beta);

gamma = 1/alpha - 1/beta;

c = zeros(size(y));

idx_left = y < gamma;

c(idx_left) = C .* exp(alpha .* (y(idx_left) - gamma)) ./ alpha^2 - y(idx_left);

c(~idx_left) = C .* exp(-beta .* (y(~idx_left) - gamma)) ./ beta^2;
end