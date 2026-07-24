function sum_price_squared = obj_fun_AB(otmSurface, eta, kappa, M, value, type)

if kappa <= 0
    sum_price_squared = 1e20;
    return;
end

p_plus = eta + sqrt((eta^2) + 1/kappa);
a = -0.5 * p_plus;

% I0 = sqrt(2*pi) * c(0; eta,kappa)
I_0 = sqrt(2*pi) * price_Bachelier_Lewis_FFT(0, 1, 1, a, 1, eta, kappa, M, value, type);

if I_0 <= 0 || ~isfinite(I_0)
    sum_price_squared = 1e20;
    return;
end

c_obs = price_Bachelier_Lewis_FFT(otmSurface.chi .* I_0, 1, 1, a, 1, eta, kappa, M, value, type);

C_model = otmSurface.DiscountFactor .* otmSurface.sigmaATM .* sqrt(otmSurface.TTM) .* (c_obs ./ I_0);

residuals = otmSurface.allCall - C_model;

sum_price_squared = sum(residuals.^2);

end