function vega = compute_vega(exotic_name, maturity, simIncrements, bump, a_norm, paramsAB, M, valueSim, typeSim, K1, num_grid_points, F_curve_all)
    
exotic_name = validatestring(exotic_name, {'CoC', 'Chooser','PoP','PV_6M','PV_12M'});
maturity = validatestring(maturity, {'6M', '12M'});

% Base AB volatilities.
sigma_base_6M = simIncrements.AB.sigma.T1;
sigma_base_12M = simIncrements.AB.sigma.T2;

% Volatility bump scenarios for the requested maturity.
switch maturity
    case '6M'
        sigma_UP_6M = sigma_base_6M + bump;
        sigma_UP_12M = sigma_base_12M;

        sigma_DW_6M = sigma_base_6M - bump;
        sigma_DW_12M = sigma_base_12M;
    case '12M'
        sigma_UP_6M = sigma_base_6M;
        sigma_UP_12M = sigma_base_12M + bump;

        sigma_DW_6M = sigma_base_6M;
        sigma_DW_12M = sigma_base_12M - bump;
end

% Simulate increments for the bumped vol scenarios.
sim_UP = build_fft_model_increments('AB', sigma_UP_6M, sigma_UP_12M,...
    a_norm, paramsAB, simIncrements.common.tGrid_d, M, valueSim, typeSim,...
    simIncrements.common.u_T0_T1, simIncrements.common.u_T1_T2, true);

sim_DW = build_fft_model_increments('AB', sigma_DW_6M, sigma_DW_12M,...
    a_norm, paramsAB, simIncrements.common.tGrid_d, M, valueSim, typeSim,...
    simIncrements.common.u_T0_T1, simIncrements.common.u_T1_T2, true);

% Dati comuni per il pricing
B_T1 = simIncrements.common.B_T1;
B_T2 = simIncrements.common.B_T2;
F_T0_T2 = simIncrements.common.F_T0_T2;
K2 = F_T0_T2; % Dal testo dell'assignment

% Chooser has no K1 strike.
if strcmp(exotic_name, 'Chooser')
    K1 = [];
end

if strcmp(exotic_name, 'PV_6M')
    F_T0_T1 = F_curve_all(3); 
    
    % S_T1 = Forward + X_T1
    S_T1_UP = F_T0_T1 + sim_UP.increments.T0_T1;
    S_T1_DW = F_T0_T1 + sim_DW.increments.T0_T1;
    
    % Price with bumped vol; K1 is the vanilla strike.
    price_UP = B_T1 * mean(max(S_T1_UP - K1, 0));
    price_DW = B_T1 * mean(max(S_T1_DW - K1, 0));

elseif strcmp(exotic_name, 'PV_12M')
    % S_T2 = Forward + X_T1 + Delta_X_T1_T2
    S_T2_UP = F_T0_T2 + sim_UP.increments.T0_T1 + sim_UP.increments.T1_T2;
    S_T2_DW = F_T0_T2 + sim_DW.increments.T0_T1 + sim_DW.increments.T1_T2;

    % Price with bumped vol; K1 is the vanilla strike.
    price_UP = B_T2 * mean(max(S_T2_UP - K1, 0));
    price_DW = B_T2 * mean(max(S_T2_DW - K1, 0));
else
    % Pricing Scenario UP
    [price_UP, ~] = price_exotic_no_grid_from_increments( ...
        exotic_name, sim_UP.increments.T0_T1, sim_UP.increments.T1_T2, ...
        B_T1, B_T2, F_T0_T2, K2, K1, num_grid_points);
    
    % Pricing Scenario DOWN
    [price_DW, ~] = price_exotic_no_grid_from_increments( ...
        exotic_name, sim_DW.increments.T0_T1, sim_DW.increments.T1_T2, ...
        B_T1, B_T2, F_T0_T2, K2, K1, num_grid_points);
end
% Central-difference vega.
vega = (price_UP - price_DW) / (2 * bump);
end