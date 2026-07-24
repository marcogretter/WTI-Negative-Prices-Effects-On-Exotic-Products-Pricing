function results = compare_MA_increment_sampling(nPaths, seed, model_MA, ...
    sigma_s_MA, sigma_t_MA, a_MA_sim, paramsMA, valueSim, typeSim, M, tGrid)
%COMPARE_MA_INCREMENT_SAMPLING
% Compares MA increment simulation methods:
%
%   1. FFT inverse CDF with estimated exponential tails
%   2. FFT inverse CDF with theoretical asymptotic tails
%   3. Analytical inverse CDF
%
% Important:
%   The FFT reconstructs the CDF of the absolutely continuous component.
%   Therefore, for a fair pathwise comparison, the analytical sampler is
%   also written as:
%       atom simulation + conditional continuous inverse CDF.
%
% Inputs:
%   nPaths       number of simulated paths
%   seed         random seed
%   model_MA     model label, usually 'MA'
%   sigma_s_MA   sigma at T1
%   sigma_t_MA   sigma at T2
%   a_MA_sim     damping parameter for FFT CDF reconstruction
%   paramsMA     struct with fields alpha, beta
%   valueSim     FFT grid parameter
%   typeSim      FFT grid type
%   M            FFT grid size parameter
%   tGrid        [T1, T2]
%
% Output:
%   results      struct with moments, timing, samples and pathwise diagnostics

    if numel(tGrid) ~= 2
        error('tGrid must be [T1, T2].');
    end

    alpha = paramsMA.alpha;
    beta  = paramsMA.beta;

    gamma = 1 / alpha - 1 / beta;

    q1 = sigma_s_MA * sqrt(tGrid(1));
    q2 = sigma_t_MA * sqrt(tGrid(2));

    if q1 <= 0 || q2 <= 0 || q2 <= q1
        error('Invalid MA scales: require 0 < q1 < q2.');
    end

    mu1 = gamma * q1;
    mu2 = gamma * q2;

    p1m = alpha / q1;
    p1p = beta  / q1;

    p2m = alpha / q2;
    p2p = beta  / q2;

    atomValue = gamma * (q2 - q1);

    % Atom and continuous probabilities
    pAtom = (p2m / p1m) * (p2p / p1p);
    pCont = 1 - pAtom;

    if pCont <= 0 || pCont >= 1
        error('Invalid continuous probability: pCont = %.12g.', pCont);
    end

    % Common random numbers:
    % U_jump decides atom vs continuous.
    % U_cont is used only conditional on being in the continuous component.
    rng(seed);

    U_jump = rand(nPaths, 1);
    U_cont = rand(nPaths, 1);

    idxCont = (U_jump < pCont);
    idxAtom = ~idxCont;

    % FFT CDF construction for the continuous component
    tic;
    cfForFFT = @(u) get_phiFFT_increment(model_MA, u, tGrid(1), tGrid(2), ...
        sigma_s_MA, sigma_t_MA, paramsMA);

    [x_grid_MA, P_grid_MA] = cdf_Lewis_FFT_from_cf( ...
        cfForFFT, a_MA_sim, M, valueSim, typeSim);
    time_fft_build = toc;

    % Inverse CDF with estimated tails
    tic;
    invCDF_est = build_inverse_cdf_spline(x_grid_MA, P_grid_MA);
    time_spline_est = toc;

    tic;
    X_est = atomValue * ones(nPaths, 1);
    X_est(idxCont) = simulated_increments(U_cont(idxCont), invCDF_est);
    time_est_sampling = toc;

    % Inverse CDF with theoretical tails
    tic;
    invCDF_theory = build_inverse_cdf_spline(x_grid_MA, P_grid_MA, p2m, p2p);
    time_spline_theory = toc;

    tic;
    X_theory = atomValue * ones(nPaths, 1);
    X_theory(idxCont) = simulated_increments(U_cont(idxCont), invCDF_theory);
    time_theory_sampling = toc;

    % Analytical inverse CDF, written with the same atom/continuous split
    tic;
    X_analytic = atomValue * ones(nPaths, 1);
    X_analytic(idxCont) = sample_increment_MA_analytic_continuous( ...
        U_cont(idxCont), p1m, p1p, p2m, p2p, mu1, mu2);
    time_analytic = toc;

    % Moments
    mom_est      = first_four_moments(X_est);
    mom_theory   = first_four_moments(X_theory);
    mom_analytic = first_four_moments(X_analytic);

    % Pathwise diagnostics
    diff_est    = X_est    - X_analytic;
    diff_theory = X_theory - X_analytic;

    stats_est    = pathwise_abs_stats(diff_est);
    stats_theory = pathwise_abs_stats(diff_theory);

    % Store results
    results.nPaths = nPaths;
    results.seed   = seed;

    results.params.alpha = alpha;
    results.params.beta  = beta;
    results.params.gamma = gamma;

    results.params.q1 = q1;
    results.params.q2 = q2;

    results.params.mu1 = mu1;
    results.params.mu2 = mu2;

    results.params.p1m = p1m;
    results.params.p1p = p1p;
    results.params.p2m = p2m;
    results.params.p2p = p2p;

    results.pCont = pCont;
    results.pAtom = pAtom;
    results.empiricalPCont = mean(idxCont);
    results.empiricalPAtom = mean(idxAtom);
    results.atomValue = atomValue;

    results.lambda.estimatedMinus = invCDF_est.lambdaMinus;
    results.lambda.estimatedPlus  = invCDF_est.lambdaPlus;
    results.lambda.theoryMinus    = invCDF_theory.lambdaMinus;
    results.lambda.theoryPlus     = invCDF_theory.lambdaPlus;

    results.samples.estimatedTails   = X_est;
    results.samples.theoreticalTails = X_theory;
    results.samples.analyticCDF      = X_analytic;

    results.time.estimatedTails   = time_fft_build + time_spline_est    + time_est_sampling;
    results.time.theoreticalTails = time_fft_build + time_spline_theory + time_theory_sampling;
    results.time.analyticCDF      = time_analytic;

    results.speedup.analyticVsEstimated   = results.time.estimatedTails   / time_analytic;
    results.speedup.analyticVsTheoretical = results.time.theoreticalTails / time_analytic;

    results.moments.names = {'Mean'; 'Variance'; 'Skewness'; 'Kurtosis'};
    results.moments.estimatedTails   = mom_est;
    results.moments.theoreticalTails = mom_theory;
    results.moments.analyticCDF      = mom_analytic;

    results.moments.absDiffEstimatedVsAnalytic = abs(mom_est - mom_analytic);
    results.moments.absDiffTheoryVsAnalytic    = abs(mom_theory - mom_analytic);

    results.pathwise.estimatedVsAnalytic   = stats_est;
    results.pathwise.theoreticalVsAnalytic = stats_theory;

    fprintf('\nMA increment pathwise comparison vs analytic CDF:\n');
    fprintf('Estimated tails   | mean abs diff = %.6e | q99 = %.6e | q999 = %.6e | max = %.6e\n', ...
        stats_est.meanAbsDiff, stats_est.q99AbsDiff, stats_est.q999AbsDiff, stats_est.maxAbsDiff);
    fprintf('Theoretical tails | mean abs diff = %.6e | q99 = %.6e | q999 = %.6e | max = %.6e\n\n', ...
        stats_theory.meanAbsDiff, stats_theory.q99AbsDiff, stats_theory.q999AbsDiff, stats_theory.maxAbsDiff);

end


function phiFFT = get_phiFFT_increment(model_MA, u, s, t, sigma_s, sigma_t, paramsMA)

    inc = cf_increment_model(model_MA, u, s, t, sigma_s, sigma_t, paramsMA);
    phiFFT = inc.phiFFT;

end


function X = sample_increment_MA_analytic_continuous(U, p1m, p1p, p2m, p2p, mu1, mu2)
% Samples the continuous component of the MA increment X_t - X_s,
% conditional on not being in the atom.

    U = U(:);
    U = min(max(U, realmin), 1 - eps);

    m = mu2 - mu1;

    A = (1 - p2m / p1m) * ...
        (p2p / p1p + ...
        (1 - p2p / p1p) * p2p / (p2m + p2p));

    Btail = (1 - p2p / p1p) * ...
        (p2m / p1m + ...
        (1 - p2m / p1m) * p2m / (p2m + p2p));

    pAtom = (p2m / p1m) * (p2p / p1p);
    pCont = 1 - pAtom;

    leftMassCond = A / pCont;

    X = zeros(size(U));

    idxLeft  = U < leftMassCond;
    idxRight = ~idxLeft;

    % Left branch:
    % U = A * exp(p2m * (x - m)) / pCont
    X(idxLeft) = m + log((U(idxLeft) .* pCont) ./ A) ./ p2m;

    % Right branch:
    % U = 1 - Btail * exp(-p2p * (x - m)) / pCont
    X(idxRight) = m - log((pCont .* (1 - U(idxRight))) ./ Btail) ./ p2p;

end


function m = first_four_moments(X)

    X = X(:);

    m = [
        mean(X);
        var(X, 1);
        skewness(X, 1);
        kurtosis(X, 1)
    ];

end


function s = pathwise_abs_stats(diffX)

    absDiff = abs(diffX(:));

    s.meanAbsDiff = mean(absDiff);
    s.medianAbsDiff = median(absDiff);
    s.q99AbsDiff = quantile(absDiff, 0.99);
    s.q999AbsDiff = quantile(absDiff, 0.999);
    s.maxAbsDiff = max(absDiff);

end