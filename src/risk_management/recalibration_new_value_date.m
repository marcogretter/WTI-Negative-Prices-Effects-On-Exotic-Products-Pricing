function [prices] = recalibration_new_value_date(projectRoot,valuationDate,K2,point4Params,resetDates)
% Extract data on new value date
mktData_new = extract_data(projectRoot,valuationDate);

% Save the plots done until now
figsBefore = findall(0, 'Type', 'figure');
% Starting calibration on new value date
curves_new = point1_calibrate_market_curves(mktData_new.snapshot, mktData_new.futureExpiries, mktData_new.valuationDate);
figsAfter = findall(0, 'Type', 'figure');
figsToClose = setdiff(figsAfter, figsBefore);
close(figsToClose);

fftParams.M     = 15;
fftParams.value = -300;
fftParams.type  = 'x1';
vol_new = point2_calibrate_vol_surface(mktData_new.snapshot, mktData_new.futureExpiries, mktData_new.valuationDate, curves_new, fftParams);

calibDiag_new = compute_point2_calibration_diagnostics(vol_new, fftParams);

% Defining parameters for simulation:
simParams.M     = 15;
simParams.valueSim = -300;
simParams.typeSim  = 'x1';
simParams.N_sim =5000000;
% Residual maturities from the new valuation date to the original T1,T2.
tGrid = yearfrac(valuationDate, resetDates, 3);

if any(tGrid <= 0)
    error('At least one reset date is not after the new valuation date.');
end

% Interpolate the newly calibrated ATM volatilities on the residual
% maturities. This replaces run_point3c_simulate_increments.
sigmaATMByMaturity_new = vol_new.sigmaATMByMaturity(:);
TTM_new = curves_new.TTM(:);
all_TTM_new = [NaN;TTM_new];

validATM = isfinite(sigmaATMByMaturity_new) & sigmaATMByMaturity_new > 0 & isfinite(all_TTM_new) & all_TTM_new > 0;

sigmaATMGrid = interp1(all_TTM_new(validATM),sigmaATMByMaturity_new(validATM),tGrid(:),'linear','extrap');

% Keep the same structure expected by build_simIncrements.
sim3c_est = struct();
sim3c_est.tGrid = tGrid(:);
sim3c_est.sigmaATMGrid = sigmaATMGrid(:);

simIncrements_new = build_simIncrements(mktData_new.snapshot,curves_new,sim3c_est.tGrid,sim3c_est.sigmaATMGrid,vol_new,simParams);

% Pricing exotics at new value date:
B_T1_new    = simIncrements_new.common.B_T1;
B_T2_new    = simIncrements_new.common.B_T2;
F_T0_T2_new = simIncrements_new.common.F_T0_T2;

[prices.price_exotic_CoC,~] = price_exotic_no_grid_from_increments( ...
    'CoC', simIncrements_new.AB.increments.T0_T1, simIncrements_new.AB.increments.T1_T2, ...
    B_T1_new, B_T2_new, F_T0_T2_new, K2, point4Params.K1_CoC, point4Params.num_grid_points);

[prices.price_exotic_PoP,~] = price_exotic_no_grid_from_increments( ...
    'PoP', simIncrements_new.AB.increments.T0_T1, simIncrements_new.AB.increments.T1_T2, ...
    B_T1_new, B_T2_new, F_T0_T2_new, K2, point4Params.K1_PoP, point4Params.num_grid_points);

[prices.price_exotic_Chooser,~] = price_exotic_no_grid_from_increments( ...
    'Chooser', simIncrements_new.AB.increments.T0_T1, simIncrements_new.AB.increments.T1_T2, ...
    B_T1_new, B_T2_new, F_T0_T2_new, K2, [], point4Params.num_grid_points);
% Outputs useful for Risk Management P&L
prices.valuationDate = valuationDate;
prices.mktData_new   = mktData_new;
prices.curves_new    = curves_new;
prices.vol_new       = vol_new;
prices.calibDiag_new = calibDiag_new;

prices.resetDates    = resetDates;
prices.tGrid         = tGrid;
prices.sigmaATMGrid  = sigmaATMGrid;

prices.B_T1_new      = B_T1_new;
prices.B_T2_new      = B_T2_new;
prices.F_T0_T2_new   = F_T0_T2_new;

end