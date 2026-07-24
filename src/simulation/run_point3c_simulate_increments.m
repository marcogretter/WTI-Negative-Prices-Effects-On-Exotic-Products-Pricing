function sim3c = run_point3c_simulate_increments(valuationDate, sigmaATMByMaturity, TTM, vol, simParams, useTheoreticalTails)
% Simulate the T1→T2 conditional increment for AB, GL, and MA models.
% Computes X(T2)−X(T1) for each model using Lewis-FFT + inverse CDF
% sampling. Supports two tail-extrapolation modes:
%   useTheoreticalTails = false (default)
%       Tail coefficients estimated from FFT CDF grid endpoints.
%       Can underestimate lambdaPlus on a coarse grid, producing an
%       artificially heavy right tail.
%
%   useTheoreticalTails = true
%       Theoretical/asymptotic tail coefficients passed explicitly to
%       build_inverse_cdf_spline, stabilising the right tail.
% Inputs:
%   valuationDate       datetime scalar (converted to datenum internally)
%   sigmaATMByMaturity  ATM vol vector (NaN where estimation failed)
%   TTM                 year fractions for valid maturities, matching
%                       the non-NaN entries of sigmaATMByMaturity
%   vol                 struct from point2_calibrate_vol_surface
%   simParams           struct: M, valueSim, typeSim, N_sim
%   useTheoreticalTails logical (default false)
%
% Output struct sim3c:
%   useTheoreticalTails, resetDates, tGrid, sigmaATMGrid
%   AB: sigma_s, sigma_t, lambdaMinus, lambdaPlus, invCDF, increments
%   GL: sigma_s, sigma_t, lambdaMinus, lambdaPlus, invCDF, increments
%   MA: sigma_s, sigma_t, lambdaMinus, lambdaPlus, invCDF, increments,
%       pCont, atomHandling

if nargin < 6
    useTheoreticalTails = false;
end

M        = simParams.M;
valueSim = simParams.valueSim;
typeSim  = simParams.typeSim;
N_sim    = simParams.N_sim;

% Reset dates and tGrid
valuationDate_dn = datenum(valuationDate);
resetDatesRaw = [ ...
    addtodate(valuationDate_dn,  6, 'month'); ...
    addtodate(valuationDate_dn, 12, 'month')];
resetDates = busdate(resetDatesRaw, 'follow');
tGrid      = yearfrac(valuationDate_dn, resetDates, 3);

T1 = tGrid(1);
T2 = tGrid(2);

% sigmaATMGrid at reset dates
idxATM       = ~isnan(sigmaATMByMaturity);
sigmaATM_vec = sigmaATMByMaturity(idxATM);
sigmaATMGrid = interp1(TTM, sigmaATM_vec, tGrid, 'linear');

% Random draws — fixed seed for reproducibility
rng(1);
u = rand(N_sim, 1);
B_MA = rand(N_sim, 1);


% AB: T1→T2 conditional increment
sigma_s_AB = sigmaATMGrid(1) / vol.AB.I0;
sigma_t_AB = sigmaATMGrid(2) / vol.AB.I0;
% CDF analyticity strip requires 0 < a < pPlusT2. a_norm is the
% call-pricing damping shift (negative); negating it gives +0.5*pPlusT2.
a_AB       = -vol.AB.a_norm / (sigma_t_AB * sqrt(T2));

cfAB = @(u) cf_increment_model('AB', u, T1, T2, sigma_s_AB, sigma_t_AB, vol.AB.params).phiFFT;
[xg_AB, Pg_AB] = cdf_Lewis_FFT_from_cf(cfAB, a_AB, M, valueSim, typeSim);

if useTheoreticalTails
    pPlus_AB  = vol.AB.eta + sqrt(vol.AB.eta^2 + 1/vol.AB.kappa);
    pMinus_AB = - vol.AB.eta + sqrt(vol.AB.eta^2 + 1/vol.AB.kappa);
    q_t_AB    = sigma_t_AB * sqrt(T2);
    lm_AB     = pMinus_AB / q_t_AB;
    lp_AB     = pPlus_AB  / q_t_AB;
    invCDF_AB = build_inverse_cdf_spline(xg_AB, Pg_AB, lm_AB, lp_AB);
else
    invCDF_AB = build_inverse_cdf_spline(xg_AB, Pg_AB);
    lm_AB     = invCDF_AB.lambdaMinus;
    lp_AB     = invCDF_AB.lambdaPlus;
end

increments_AB = simulated_increments(u, invCDF_AB);

% GL: T1→T2 conditional increment
sigma_s_GL = sigmaATMGrid(1) / vol.GL.I0;
sigma_t_GL = sigmaATMGrid(2) / vol.GL.I0;
% CDF analyticity strip requires 0 < a < pPlus. a_norm is the
% call-pricing damping shift (negative); negating it gives +0.5*pPlusT2.
a_GL       = -vol.GL.a_norm / (sigma_t_GL * sqrt(T2));

cfGL = @(u) cf_increment_model('GL', u, T1, T2, sigma_s_GL, sigma_t_GL, vol.GL.params).phiFFT;
[xg_GL, Pg_GL] = cdf_Lewis_FFT_from_cf(cfGL, a_GL, M, valueSim, typeSim);

if useTheoreticalTails
    q_t_GL = sigma_t_GL * sqrt(T2);
    lm_GL  = vol.GL.alpha / q_t_GL;
    lp_GL  = vol.GL.beta  / q_t_GL;
    invCDF_GL = build_inverse_cdf_spline(xg_GL, Pg_GL, lm_GL, lp_GL);
else
    invCDF_GL = build_inverse_cdf_spline(xg_GL, Pg_GL);
    lm_GL     = invCDF_GL.lambdaMinus;
    lp_GL     = invCDF_GL.lambdaPlus;
end

increments_GL = simulated_increments(u, invCDF_GL);

% MA: T1→T2 conditional increment with atom
% manage_MA_increments is not called here because it does not return the
% invCDF object needed for the output struct. The atom logic is inlined.
sigma_s_MA = sigmaATMGrid(1) / vol.MA.I0;
sigma_t_MA = sigmaATMGrid(2) / vol.MA.I0;
% CDF analyticity strip requires 0 < a < pPlus. a_norm is the
% call-pricing damping shift (negative); negating it gives +0.5*pPlusT2.
a_MA       = -vol.MA.a_norm / (sigma_t_MA * sqrt(T2));

cfMA = @(u) cf_increment_model('MA', u, T1, T2, sigma_s_MA, sigma_t_MA, vol.MA.params).phiFFT;
[xg_MA, Pg_MA] = cdf_Lewis_FFT_from_cf(cfMA, a_MA, M, valueSim, typeSim);

if useTheoreticalTails
    q_t_MA = sigma_t_MA * sqrt(T2);
    lm_MA  = vol.MA.alpha / q_t_MA;
    lp_MA  = vol.MA.beta  / q_t_MA;
    invCDF_MA = build_inverse_cdf_spline(xg_MA, Pg_MA, lm_MA, lp_MA);
else
    invCDF_MA = build_inverse_cdf_spline(xg_MA, Pg_MA);
    lm_MA     = invCDF_MA.lambdaMinus;
    lp_MA     = invCDF_MA.lambdaPlus;
end

p_cont_MA    = 1 - (sigma_s_MA^2 * T1) / (sigma_t_MA^2 * T2);
gammaLoc_MA  = 1/vol.MA.alpha - 1/vol.MA.beta;
atomValue_MA = gammaLoc_MA * (sigma_t_MA * sqrt(T2) - sigma_s_MA * sqrt(T1));

Bern_MA       = double(B_MA < p_cont_MA);
increments_MA = atomValue_MA * ones(N_sim, 1);
idxCont       = (Bern_MA == 1);
if any(idxCont)
    increments_MA(idxCont) = simulated_increments(u(idxCont), invCDF_MA);
end

% Output struct
sim3c.useTheoreticalTails = useTheoreticalTails;
sim3c.resetDates          = resetDates;
sim3c.tGrid               = tGrid;
sim3c.sigmaATMGrid        = sigmaATMGrid;

sim3c.AB.sigma_s     = sigma_s_AB;
sim3c.AB.sigma_t     = sigma_t_AB;
sim3c.AB.lambdaMinus = lm_AB;
sim3c.AB.lambdaPlus  = lp_AB;
sim3c.AB.invCDF      = invCDF_AB;
sim3c.AB.increments  = increments_AB;

sim3c.GL.sigma_s     = sigma_s_GL;
sim3c.GL.sigma_t     = sigma_t_GL;
sim3c.GL.lambdaMinus = lm_GL;
sim3c.GL.lambdaPlus  = lp_GL;
sim3c.GL.invCDF      = invCDF_GL;
sim3c.GL.increments  = increments_GL;

sim3c.MA.sigma_s                = sigma_s_MA;
sim3c.MA.sigma_t                = sigma_t_MA;
sim3c.MA.lambdaMinus            = lm_MA;
sim3c.MA.lambdaPlus             = lp_MA;
sim3c.MA.invCDF                 = invCDF_MA;
sim3c.MA.increments             = increments_MA;
sim3c.MA.pCont                  = p_cont_MA;
sim3c.MA.atomHandling.atomValue = atomValue_MA;
sim3c.MA.atomHandling.gammaLoc  = gammaLoc_MA;
sim3c.MA.atomHandling.Bernoulli = Bern_MA;

end
