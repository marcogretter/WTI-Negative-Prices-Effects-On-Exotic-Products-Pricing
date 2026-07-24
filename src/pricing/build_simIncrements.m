function simIncrements = build_simIncrements(snapshot, curves, tGrid, sigmaATMGrid, vol, simParams)
% Build simulated increments for AB, GL and MA models.
%
% Inputs:
%   snapshot      struct array — snapshot(1).valuationDate used for T0
%   curves        struct from point1_calibrate_market_curves
%   tGrid         2×1 year fractions [T1; T2] from run_point3c_simulate_increments
%   sigmaATMGrid  2×1 ATM vols at [T1; T2]
%   vol           struct from point2_calibrate_vol_surface
%   simParams     struct: M, valueSim, typeSim

M        = simParams.M;
valueSim = simParams.valueSim;
typeSim  = simParams.typeSim;

futureExpiries_curve = curves.futureExpiries_curve;
B_curve              = curves.B_curve;
F_curve              = curves.F_curve;
TTM                  = curves.TTM;

I0_AB    = vol.AB.I0;      aAB_norm = vol.AB.a_norm;   paramsAB = vol.AB.params;
I0_GL    = vol.GL.I0;      a_GL_norm = vol.GL.a_norm;  paramsGL = vol.GL.params;
I0_MA    = vol.MA.I0;      a_MA_norm = vol.MA.a_norm;  paramsMA = vol.MA.params;

% Dates and market quantities

T0 = snapshot(1).valuationDate;
T1 = dateAddMonth(T0, 6);
T2 = dateAddMonth(T0, 12);

B_T1 = get_discount_factor_by_zero_rates_linear_interp(T0, T1, futureExpiries_curve, B_curve);

B_T2 = get_discount_factor_by_zero_rates_linear_interp(T0, T2, futureExpiries_curve, B_curve);

tGrid_d = [0; tGrid];

F_T0_T2 = interp1(TTM, F_curve, tGrid_d(3));

% Simulation setup
if isfield(simParams, 'N_sim')
    N_sim = simParams.N_sim;
else
    N_sim = 1000000';
end

if isfield(simParams, 'seed')
    seed = simParams.seed;
else
    seed = 1;
end

if isfield(simParams, 'useTheoreticalTails')
    useTheoreticalTails = simParams.useTheoreticalTails;
else
    useTheoreticalTails = true;
end

rng(seed);

u_T0_T1 = rand(N_sim, 1);
u_T1_T2 = rand(N_sim, 1);

% Store common quantities

simIncrements = struct();

simIncrements.common.T0 = T0;
simIncrements.common.T1 = T1;
simIncrements.common.T2 = T2;

simIncrements.common.B_T1    = B_T1;
simIncrements.common.B_T2    = B_T2;
simIncrements.common.B12     = B_T2 / B_T1;
simIncrements.common.F_T0_T2 = F_T0_T2;
simIncrements.common.tGrid_d = tGrid_d;

simIncrements.common.N_sim    = N_sim;
simIncrements.common.seed     = 1;
simIncrements.common.u_T0_T1  = u_T0_T1;
simIncrements.common.u_T1_T2  = u_T1_T2;

% AB increments

sigma_T1_AB = sigmaATMGrid(1) / I0_AB;
sigma_T2_AB = sigmaATMGrid(2) / I0_AB;

simIncrements.AB = build_fft_model_increments('AB', sigma_T1_AB, sigma_T2_AB, aAB_norm, paramsAB, tGrid_d, M, valueSim, typeSim,u_T0_T1, u_T1_T2, useTheoreticalTails);

% GL increments

sigma_T1_GL = sigmaATMGrid(1) / I0_GL;
sigma_T2_GL = sigmaATMGrid(2) / I0_GL;

simIncrements.GL = build_fft_model_increments('GL', sigma_T1_GL, sigma_T2_GL, a_GL_norm, paramsGL, tGrid_d, M, valueSim, typeSim, u_T0_T1, u_T1_T2, useTheoreticalTails);

% MA increments

sigma_T1_MA = sigmaATMGrid(1) / I0_MA;
sigma_T2_MA = sigmaATMGrid(2) / I0_MA;

simIncrements.MA = build_MA_model_increments(sigma_T1_MA, sigma_T2_MA, a_MA_norm,paramsMA, tGrid_d, M, valueSim, typeSim, u_T0_T1, u_T1_T2, useTheoreticalTails);

% Store model-specific calibrated sigmas

simIncrements.AB.sigma.T1 = sigma_T1_AB;
simIncrements.AB.sigma.T2 = sigma_T2_AB;

simIncrements.GL.sigma.T1 = sigma_T1_GL;
simIncrements.GL.sigma.T2 = sigma_T2_GL;

simIncrements.MA.sigma.T1 = sigma_T1_MA;
simIncrements.MA.sigma.T2 = sigma_T2_MA;

simIncrements.common.seed = seed;
simIncrements.common.useTheoreticalTails = useTheoreticalTails;
fprintf('\nSimulated increments created.\n');

end
