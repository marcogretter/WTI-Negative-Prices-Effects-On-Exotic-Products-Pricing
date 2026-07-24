function curves = point1_calibrate_market_curves(snapshot, futureExpiries, valuationDate)
%futureExpiries  datetime vector of all maturity dates
%   valuationDate   datetime scalar
%
% Output struct curves fields:
%   B_curve_all, F_curve_all, R2_curve_all   all maturities (NaN where invalid)
%   B_curve, F_curve, R Point 1: calibrate discount factors, forwards, zero rates, and absolute dividends.
%
% Inputs:
%   snapshot        struct array with fields K, C, P per maturity
%   2_curve, TTM          valid maturities only
%   zero_rates                               continuously compounded, valid maturities
%   futureExpiries_all                       datetime N×1, all expiries
%   futureExpiries_curve                     datetime M×1, valid maturities only

% Regression to find DFs and Fs (point 1 of the Assignment).
[B_curve_all, F_curve_all, R2_curve_all] = calibrate_discount_curve( ...
    {snapshot.K}, {snapshot.C}, {snapshot.K}, {snapshot.P});

futureExpiries_all = futureExpiries(:);
TTM_all = yearfrac(valuationDate, futureExpiries_all, 3);   % ACT/365

% The first expiry date of the future file already expired at the valuation
% Date
validMaturities = TTM_all > 0 & isfinite(B_curve_all) & B_curve_all > 0 & isfinite(F_curve_all);

B_curve = B_curve_all(validMaturities);
F_curve = F_curve_all(validMaturities);
R2_curve = R2_curve_all(validMaturities);
futureExpiries_curve = futureExpiries_all(validMaturities);
TTM = TTM_all(validMaturities);

zero_rates = plot_calibration_df_and_forwards(B_curve, F_curve, futureExpiries_curve, valuationDate);

% Absolute dividends / implied carry between consecutive expiries
% Formula:
% D_{j-1,j}^{T_j} = (B(0,T_{j-1}) / B(0,T_j)) * F(0,T_{j-1}) - F(0,T_j)
% synthetic forwards from put-call parity.

capitalization_factor = B_curve(1:end-1) ./ B_curve(2:end);
absolute_dividend = capitalization_factor .* F_curve(1:end-1) - F_curve(2:end);

delta_T = diff(TTM);
absolute_dividend_intensity = absolute_dividend ./ delta_T;
figure('Name', 'Absolute Dividends / Implied Carry', 'NumberTitle', 'off', 'Color', 'w');

subplot(2,1,1);
plot(TTM(2:end), absolute_dividend, '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
xlabel('Final Time to Maturity (years)');
ylabel('Absolute Dividend');
title('Implied absolute dividends between consecutive maturities');

subplot(2,1,2);
stairs(TTM(2:end), absolute_dividend_intensity, 'LineWidth', 1.5);
hold on;
plot(TTM(2:end), absolute_dividend_intensity, 's', 'MarkerSize', 8);
grid on;
xlabel('Final Time to Maturity (years)');
ylabel('Absolute Dividend Rate');
title('Piecewise-constant absolute dividend rates');

curves.B_curve_all          = B_curve_all;
curves.F_curve_all          = F_curve_all;
curves.R2_curve_all         = R2_curve_all;
curves.B_curve              = B_curve;
curves.F_curve              = F_curve;
curves.R2_curve             = R2_curve;
curves.TTM                  = TTM;
curves.zero_rates           = zero_rates;
curves.futureExpiries_all   = futureExpiries_all;
curves.futureExpiries_curve = futureExpiries_curve;
end