function [alphaGL, betaGL] = calibrate_alpha_beta_GL(otmSurface, alpha_0, beta_0, M, value, type)
% Calibrate Generalized Logistic parameters alpha and beta
% using the same cascade objective structure as AB.
%
% Requires:
%   obj_fun_GL.m
%   price_Bachelier_Lewis_FFT.m modified with model flag 'GL'
%   generalizedLogisticCf.m
%   complexLogGammaLanczos.m

p0 = [alpha_0, beta_0];

options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000, 'TolX', 1e-6, 'TolFun', 1e-8);

p_results = fminsearch(@(p) obj_fun_GL(otmSurface, p(1), p(2), M, value, type), p0, options);

alphaGL = p_results(1);
betaGL  = p_results(2);
end