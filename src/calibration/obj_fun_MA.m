function sum_price_squared = obj_fun_MA(otmSurface, alpha, beta)

if alpha <= 0 || beta <= 0 || ~isfinite(alpha) || ~isfinite(beta)
    sum_price_squared = 1e20;
    return;
end

% I0 = sqrt(2*pi) * c_MA(0)
c0 = call_MA_normalized(0, alpha, beta);

I_0 = sqrt(2*pi) * c0;

if I_0 <= 0 || ~isfinite(I_0)
    sum_price_squared = 1e20;
    return;
end

% y_i = I0 * chi_i
y_obs = I_0 .* otmSurface.chi;

% c_MA(I0 * chi_i)
c_obs = call_MA_normalized(y_obs, alpha, beta);

if any(~isfinite(c_obs))
    sum_price_squared = 1e20;
    return;
end

% C_model = B * sigmaATM * sqrt(T) * c(I0 chi) / I0
C_model = otmSurface.DiscountFactor .* otmSurface.sigmaATM .* sqrt(otmSurface.TTM) .* (c_obs ./ I_0);

residuals = otmSurface.allCall - C_model;

sum_price_squared = sum(residuals.^2);

end