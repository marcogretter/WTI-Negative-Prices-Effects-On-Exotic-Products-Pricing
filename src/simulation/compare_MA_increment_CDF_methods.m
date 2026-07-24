function cmp = compare_MA_increment_CDF_methods(vol, sim3c, simParams)
% Compare MA increment simulation methods: FFT-estimated tails,
% FFT-theoretical tails, and closed-form analytical inverse CDF.
%
% Inputs:
%   vol       struct from point2_calibrate_vol_surface
%   sim3c     struct from run_point3c_simulate_increments
%   simParams struct: M, valueSim, typeSim, N_sim
%
% Delegates computation to compare_MA_increment_sampling.
% Returns a struct with per-method increments, timing, moments table,
% and absolute-difference table vs the analytic reference.

SEED = 2;   % dedicated seed, independent of main simulation (seed 1)

sigma_s_MA = sim3c.MA.sigma_s;
sigma_t_MA = sim3c.MA.sigma_t;
tGrid      = sim3c.tGrid;

T2       = tGrid(2);
% CDF analyticity strip requires 0 < a < lambdaPlus. a_norm is the
% call-pricing damping shift (negative); negating it gives +0.5*lambdaPlus.
a_MA_sim = -vol.MA.a_norm / (sigma_t_MA * sqrt(T2));

raw = compare_MA_increment_sampling(simParams.N_sim, SEED, 'MA',sigma_s_MA, sigma_t_MA, a_MA_sim, vol.MA.params, ...
    simParams.valueSim, simParams.typeSim, simParams.M, tGrid);

% Per-method fields
cmp.analytic.increments   = raw.samples.analyticCDF;
cmp.analytic.time         = raw.time.analyticCDF;

cmp.fft_est.increments    = raw.samples.estimatedTails;
cmp.fft_est.time          = raw.time.estimatedTails;

cmp.fft_theory.increments = raw.samples.theoreticalTails;
cmp.fft_theory.time       = raw.time.theoreticalTails;

% Moments table
methodNames = {'FFT-Estimated'; 'FFT-Theoretical'; 'Analytic-CDF'};
times = [ ...
    raw.time.estimatedTails; ...
    raw.time.theoreticalTails; ...
    raw.time.analyticCDF];

costVsAnalytic = times ./ raw.time.analyticCDF;

mom = [ ...
    raw.moments.estimatedTails, ...
    raw.moments.theoreticalTails, ...
    raw.moments.analyticCDF]';   % 3×4

cmp.momentsTable = table( ...
    methodNames, times, costVsAnalytic, ...
    mom(:,1), mom(:,2), mom(:,3), mom(:,4), ...
    'VariableNames', { ...
        'Method', 'TotalTimeSeconds', 'CostVsAnalytic', ...
        'Mean', 'Variance', 'Skewness', 'Kurtosis' ...
    } ...
);
% Absolute difference table (analytic as reference)
metricNames = {'Mean'; 'Variance'; 'Skewness'; 'Kurtosis'};

cmp.diffTable = table( ...
    metricNames, ...
    raw.moments.absDiffEstimatedVsAnalytic, ...
    raw.moments.absDiffTheoryVsAnalytic, ...
    'VariableNames', { ...
        'Metric', ...
        'AbsDiff_FFTEstimated_vs_Analytic', ...
        'AbsDiff_FFTTheoretical_vs_Analytic' ...
    } ...
);
% Settings
cmp.settings.N_sim    = simParams.N_sim;
cmp.settings.seed     = SEED;
cmp.settings.T1       = tGrid(1);
cmp.settings.T2       = T2;
cmp.settings.sigma_s  = sigma_s_MA;
cmp.settings.sigma_t  = sigma_t_MA;
cmp.settings.alpha    = vol.MA.alpha;
cmp.settings.beta     = vol.MA.beta;
cmp.settings.M        = simParams.M;
cmp.settings.valueSim = simParams.valueSim;
cmp.settings.typeSim  = simParams.typeSim;

end
