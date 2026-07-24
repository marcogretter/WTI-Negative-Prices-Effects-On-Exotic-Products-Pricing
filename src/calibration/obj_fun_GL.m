function sum_price_squared = obj_fun_GL(otmSurface, alphaGL, betaGL, M, value, type)

if alphaGL <= 0 || betaGL <= 0 || ~isfinite(alphaGL) || ~isfinite(betaGL)
    sum_price_squared = 1e20;
    return;
end

% GL analyticity strip: a in (-betaGL, alphaGL).
% Choose negative shift, so R_a = 0 in Lewis.
a = -0.5 * betaGL;

% c_GL(0)
c0 = price_Bachelier_Lewis_FFT(0, 1, 1, a, 1, alphaGL, betaGL, M, value, type, 'GL');

I_0 = sqrt(2*pi) * c0;

if I_0 <= 0 || ~isfinite(I_0)
    sum_price_squared = 1e20;
    return;
end

% y_i = I0 * chi_i
y_obs = I_0 .* otmSurface.chi;

% c_GL(I0 * chi_i)
c_obs = price_Bachelier_Lewis_FFT(y_obs, 1, 1, a, 1, alphaGL, betaGL, M, value, type, 'GL');

if any(~isfinite(c_obs))
    sum_price_squared = 1e20;
    return;
end

C_model = otmSurface.DiscountFactor .* otmSurface.sigmaATM .* sqrt(otmSurface.TTM) .* (c_obs ./ I_0);

residuals = otmSurface.allCall - C_model;

sum_price_squared = sum(residuals.^2);
end