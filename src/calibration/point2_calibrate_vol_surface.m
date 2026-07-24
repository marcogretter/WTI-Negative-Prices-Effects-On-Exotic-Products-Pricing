function vol = point2_calibrate_vol_surface(snapshot, futureExpiries, valuationDate, curves, fftParams)
% Point 2: build the Fwd-OTM surface, estimate sigma_ATM, and calibrate
% the three models (AB, MA, GL).
%
% Inputs:
%   snapshot        struct array (one element per maturity)
%   futureExpiries  datetime N×1, all expiry dates
%   valuationDate   datetime scalar
%   curves          output struct from point1_calibrate_market_curves
%   fftParams       struct with fields M, value, type
%
% Output struct vol fields:
%   otmSurface           enriched Fwd-OTM table
%   sigmaATMByMaturity   ATM vol vector (NaN where not available)
%   AB.eta, AB.kappa, AB.a_norm, AB.I0, AB.params
%   MA.alpha, MA.beta,  MA.a_norm, MA.I0, MA.params
%   GL.alpha, GL.beta,  GL.a_norm, GL.I0, GL.params

B_curve_all = curves.B_curve_all;
F_curve_all = curves.F_curve_all;

M     = fftParams.M;
value = fftParams.value;
type  = fftParams.type;

% Build Fwd-OTM surface
otmSurface = buildFwdOtmSurfacePoint2(snapshot, futureExpiries, B_curve_all, F_curve_all);

% Estimate sigma_ATM by maturity
[sigmaATMByMaturity, otmSurface] = estimateSigmaATMByMaturity(otmSurface);

% Calibration: Additive Bachelier (AB)
eta_0   = 0.2;
kappa_0 = 1.0;
[eta, kappa] = calibrate_eta_kappa(otmSurface, eta_0, kappa_0, M, value, type);

% Calibration: Minimal Additive (MA)
maInitialGrid = [
    0.5 0.5
    0.5 1.0
    0.5 2.0
    0.5 5.0
    1.0 0.5
    1.0 1.0
    1.0 2.0
    1.0 5.0
    2.0 0.5
    2.0 1.0
    2.0 2.0
    2.0 5.0
    5.0 0.5
    5.0 1.0
    5.0 2.0
    5.0 5.0
];

alpha_0 = 0.9;
beta_0  = 1.5;

[alpha_MA, beta_MA, maStartTable] = calibrate_alpha_beta_MA(otmSurface, alpha_0, beta_0, maInitialGrid);

% Calibration: Generalized Logistic (GL)
alpha_0_GL = 0.9;
beta_0_GL  = 1.5;
[alpha_GL, beta_GL] = calibrate_alpha_beta_GL(otmSurface, alpha_0_GL, beta_0_GL, M, value, type);

% I0 and a_norm constants for each model
pPlusAB  = eta + sqrt(eta^2 + 1/kappa);
aAB_norm = -0.5 * pPlusAB;
c0_AB    = price_Bachelier_Lewis_FFT(0, 1, 1, aAB_norm, 1, eta, kappa, M, value, type, 'AB');
I0_AB    = sqrt(2*pi) * c0_AB;

a_MA_norm = -0.5 * beta_MA;
c0_MA     = call_MA_normalized(0, alpha_MA, beta_MA);
I0_MA     = sqrt(2*pi) * c0_MA;

a_GL_norm = -0.5 * beta_GL;
c0_GL     = price_Bachelier_Lewis_FFT(0, 1, 1, a_GL_norm, 1, alpha_GL, beta_GL, M, value, type, 'GL');
I0_GL     = sqrt(2*pi) * c0_GL;

% Output struct
vol.otmSurface         = otmSurface;
vol.sigmaATMByMaturity = sigmaATMByMaturity;

vol.AB.eta    = eta;
vol.AB.kappa  = kappa;
vol.AB.a_norm = aAB_norm;
vol.AB.I0     = I0_AB;
vol.AB.params.eta   = eta;
vol.AB.params.kappa = kappa;

vol.MA.alpha  = alpha_MA;
vol.MA.beta   = beta_MA;
vol.MA.a_norm = a_MA_norm;
vol.MA.I0     = I0_MA;
vol.MA.params.alpha = alpha_MA;
vol.MA.params.beta  = beta_MA;
vol.MA.multistartTable = maStartTable;

vol.GL.alpha  = alpha_GL;
vol.GL.beta   = beta_GL;
vol.GL.a_norm = a_GL_norm;
vol.GL.I0     = I0_GL;
vol.GL.params.alpha = alpha_GL;
vol.GL.params.beta  = beta_GL;

end
