function results = compare_fso_MA_MC_vs_analytic( ...
    mcPrice, mcSE, F0T2, K2, B0T2, B1T2, paramsMA, sigma_s_MA, sigma_t_MA, tGrid)

% COMPARE_FSO_MA_MC_VS_ANALYTIC
%
% Compares a pre-computed MC price of the MA forward-start option against
% the closed-form analytical price.
%
% Payoff at T2:
%
%   [ S_T2 - K2 F(T1,T2) ]^+
%
% Inputs:
%   mcPrice    scalar MC price (already discounted)
%   mcSE       MC standard error
%   F0T2       F(0,T2)
%   K2         forward-start strike multiplier
%   B0T2       B(0,T2)
%   B1T2       B(T1,T2)
%   paramsMA   struct with fields alpha, beta
%   sigma_s_MA sigma at T1
%   sigma_t_MA sigma at T2
%   tGrid      [T1, T2]

    alpha = paramsMA.alpha;
    beta  = paramsMA.beta;
    gamma = 1 / alpha - 1 / beta;

    q1  = sigma_s_MA * sqrt(tGrid(1));
    q2  = sigma_t_MA * sqrt(tGrid(2));
    mu1 = gamma * q1;
    mu2 = gamma * q2;
    p1m = alpha / q1;   p1p = beta / q1;
    p2m = alpha / q2;   p2p = beta / q2;

    priceAnalytic = price_fso_MA_analytic(F0T2, K2, B0T2, B1T2, mu1, mu2, p1m, p1p, p2m, p2p);

    absError = abs(mcPrice - priceAnalytic);
    ciLow    = mcPrice - 1.96 * mcSE;
    ciHigh   = mcPrice + 1.96 * mcSE;

    if abs(priceAnalytic) > eps
        relError = absError / priceAnalytic;
    else
        relError = NaN;
    end

    if mcSE > 0
        zScore = (mcPrice - priceAnalytic) / mcSE;
    else
        zScore = NaN;
    end

    results.K2            = K2;
    results.F0T2          = F0T2;
    results.B0T2          = B0T2;
    results.B1T2          = B1T2;
    results.priceMC       = mcPrice;
    results.seMC          = mcSE;
    results.ciLow         = ciLow;
    results.ciHigh        = ciHigh;
    results.priceAnalytic = priceAnalytic;
    results.absError      = absError;
    results.relError      = relError;
    results.zScore        = zScore;
    results.params.mu1    = mu1;   results.params.mu2 = mu2;
    results.params.p1m    = p1m;   results.params.p1p = p1p;
    results.params.p2m    = p2m;   results.params.p2p = p2p;

    results.table = table( ...
        K2, priceAnalytic, mcPrice, mcSE, ciLow, ciHigh, absError, relError, ...
        'VariableNames', {'K2','AnalyticPrice','MCPrice','MCStdError','CI95_Low','CI95_High','AbsError','RelError'} ...
    );

end
