% RunFinalProject2b
% Main of the Final Project of Group 2b of Financial Engineering.
clear all; clc; close all;

projectRoot = pwd;
cd(projectRoot);

addpath(genpath(fullfile(projectRoot, 'src')));
%% Data loading
valuationDate = datetime(2020, 6, 2);
mktData        = extract_data(projectRoot,valuationDate);
snapshot       = mktData.snapshot;
futureExpiries = mktData.futureExpiries;
valuationDate  = mktData.valuationDate;
%% Point 0: comparison of the PDFs of the models
point0 = point0_plot_pdfs();
%% Point 1: calibration discount/forward and absolute dividends
curves = point1_calibrate_market_curves(snapshot, futureExpiries, valuationDate);
%% Point 2: Calibrate the Implied Volatility Surface
fftParams.M     = 15;
fftParams.value = -300;
fftParams.type  = 'x1';
vol = point2_calibrate_vol_surface(snapshot, futureExpiries, valuationDate, curves, fftParams);
%% Point 2 diagnostics
calibDiag = compute_point2_calibration_diagnostics(vol, fftParams);
disp('Calibration quality — global:');
disp(calibDiag.global);
disp('Calibration quality — by maturity:');
disp(calibDiag.byMaturity);
disp('Calibration quality — by moneyness bucket:');
disp(calibDiag.byChiBucket);
disp('Calibration quality — by side:');
disp(calibDiag.bySide);
if isfield(calibDiag, 'maMultistartTable')
    disp('MA multi-start calibration results:');
    disp(calibDiag.maMultistartTable);
end
%% Punto 3c: T1→T2 increment simulation — estimated vs theoretical tail comparison
simParams3c.M        = fftParams.M;
simParams3c.valueSim = -300;
simParams3c.typeSim  = 'x1';
simParams3c.N_sim    = 5000000;
% Estimated tail coefficients:
sim3c_est    = run_point3c_simulate_increments(valuationDate, vol.sigmaATMByMaturity, curves.TTM, vol, simParams3c, false);
% Asymptotic tail coefficients:
sim3c_theory = run_point3c_simulate_increments(valuationDate, vol.sigmaATMByMaturity, curves.TTM, vol, simParams3c, true);

modelLabels = {'AB', 'GL', 'MA'};
tailCompTable = table( ...
    modelLabels(:), ...
    [sim3c_est.AB.lambdaMinus;  sim3c_est.GL.lambdaMinus;  sim3c_est.MA.lambdaMinus], ...
    [sim3c_theory.AB.lambdaMinus; sim3c_theory.GL.lambdaMinus; sim3c_theory.MA.lambdaMinus], ...
    [sim3c_est.AB.lambdaPlus;   sim3c_est.GL.lambdaPlus;   sim3c_est.MA.lambdaPlus], ...
    [sim3c_theory.AB.lambdaPlus;  sim3c_theory.GL.lambdaPlus;  sim3c_theory.MA.lambdaPlus], ...
    [max(sim3c_est.AB.increments);  max(sim3c_est.GL.increments);  max(sim3c_est.MA.increments)], ...
    [max(sim3c_theory.AB.increments); max(sim3c_theory.GL.increments); max(sim3c_theory.MA.increments)], ...
    [std(sim3c_est.AB.increments);  std(sim3c_est.GL.increments);  std(sim3c_est.MA.increments)], ...
    [std(sim3c_theory.AB.increments); std(sim3c_theory.GL.increments); std(sim3c_theory.MA.increments)], ...
    'VariableNames', {'Model', ...
        'lambdaMinus_Est',   'lambdaMinus_Theory', ...
        'lambdaPlus_Est',    'lambdaPlus_Theory', ...
        'maxIncrement_Est',  'maxIncrement_Theory', ...
        'stdIncrement_Est',  'stdIncrement_Theory'} ...
);
disp('Tail-mode comparison: estimated vs theoretical tail coefficients');
disp(tailCompTable);
%% Point 3c diagnostics
cdfDiag = compute_point3c_cdf_diagnostics(sim3c_est, sim3c_theory, simParams3c);
disp('Lewis-FFT CDF reconstruction quality (AB / GL / MA, estimated and theoretical):');
disp(cdfDiag.cdfQuality);
disp('Right-tail comparison: estimated vs theoretical:');
disp(cdfDiag.tailComparison);
%% Point 3e optional: MA analytical CDF comparison
maCDFComparison = compare_MA_increment_CDF_methods(vol, sim3c_est, simParams3c);
disp('MA increment — methods comparison (moments and timing):');
disp(maCDFComparison.momentsTable);
disp('Absolute differences vs analytic CDF:');
disp(maCDFComparison.diffTable);
%% Point 3d: build simulated increments for each model
simIncrements = build_simIncrements(snapshot, curves, sim3c_est.tGrid, sim3c_est.sigmaATMGrid, vol, simParams3c);

B_T1        = simIncrements.common.B_T1;
B_T2        = simIncrements.common.B_T2;
F_T0_T2     = simIncrements.common.F_T0_T2;
tGrid_d     = simIncrements.common.tGrid_d;
sigma_T1_MA = simIncrements.MA.sigma.T1;
sigma_T2_MA = simIncrements.MA.sigma.T2;
%% Point 3d: forward-start option prices
k_2      = 1;
modelNames = {'AB'; 'GL'; 'MA'};
N_sim    = simIncrements.common.N_sim;
fsoPrice = zeros(3,1);
fsoSE    = zeros(3,1);
for i = 1:numel(modelNames)
    modelName = modelNames{i};
    x_T0_T1 = simIncrements.(modelName).increments.T0_T1;
    x_T1_T2 = simIncrements.(modelName).increments.T1_T2;
    payoff = max(F_T0_T2 .* (1 - k_2) + x_T0_T1 .* (1 - k_2 .* (B_T1 ./ B_T2)) + x_T1_T2, 0);
    fsoPrice(i) = B_T2 * mean(payoff);
    fsoSE(i)    = B_T2 * std(payoff) / sqrt(N_sim);
    simIncrements.(modelName).prices.forwardStart   = fsoPrice(i);
    simIncrements.(modelName).prices.forwardStartSE = fsoSE(i);
end
disp(table(modelNames, fsoPrice, fsoSE, fsoPrice-1.96*fsoSE, fsoPrice+1.96*fsoSE, ...
    'VariableNames', {'Model','MCPrice','MCStdError','CI95_Low','CI95_High'}));
%% Point 3f (optional): MA forward-start option — MC vs analytic
fsoMA = compare_fso_MA_MC_vs_analytic(fsoPrice(3), fsoSE(3), F_T0_T2, k_2, B_T2, B_T2/B_T1, ...
    vol.MA.params, sigma_T1_MA, sigma_T2_MA, tGrid_d(2:3));

disp('MA forward-start option: MC vs analytic');
disp(fsoMA.table);
%% Point 4: exotic pricing — grid vs no-grid benchmark
point4Params.num_grid_points = 2000;
point4Params.K1_CoC          = 0.3;
point4Params.K1_PoP          = 0.3;

point4 = run_point4_exotic_pricing(simIncrements, vol, point4Params);

% Required by Risk_Management (script reads these from workspace)
ABRow            = strcmp(point4.tables.gridPrices.Model, 'AB');
initialprices.CoC_AB     = point4.tables.gridPrices{ABRow, 'CoC'};
initialprices.PoP_AB     = point4.tables.gridPrices{ABRow, 'PoP'};
initialprices.Chooser_AB = point4.tables.gridPrices{ABRow, 'Chooser'};

%% Point 5: Compare prices Closed Vs Numerical formula for MA model
disp('Comparing closed Vs Analytical CoC formula in the MA model');
disp(point4.diagnostics.MA_analytical.CoC);
disp('Comparing closed Vs Analytical PoP formula in the MA model');
disp(point4.diagnostics.MA_analytical.PoP);
disp('Comparing closed Vs Analytical Chooser formula in the MA model');
disp(point4.diagnostics.MA_analytical.Chooser);

%% Point 6: Risk Management:
point6=Risk_Management(mktData,curves,point4Params,simIncrements,vol,simParams3c,projectRoot,initialprices);
disp('Risk Management: P&L summary');
disp(point6.tables.pnlSummary);

disp('Risk Management: hedge positions');
disp(point6.tables.hedgePositions);

disp('Risk Management: Greeks summary');
disp(point6.tables.greeksSummary);

fprintf('Risk Management: Cash Flow on the 2nd June 2020: %8.4f\n',point6.CF0);