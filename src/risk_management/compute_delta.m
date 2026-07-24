function delta = compute_delta(exotic_name, simIncrements, bump, K1, num_grid_points, F_curve_all)
    
exotic_name = validatestring(exotic_name, {'CoC', 'Chooser','PoP','PV_6M','PV_12M'});

% Data extraction
B_T1 = simIncrements.common.B_T1;
B_T2 = simIncrements.common.B_T2;

% AB model increments (pre-simulated, common random numbers).
inc_T0_T1 = simIncrements.AB.increments.T0_T1;
inc_T1_T2 = simIncrements.AB.increments.T1_T2;

if strcmp(exotic_name, 'PV_6M')
    % 6M forward price.
    F_base_6M = F_curve_all(3);

    % Forward bump scenarios (up and down).
    F_UP = F_base_6M + bump;
    F_DW = F_base_6M - bump;

    % S_T1 = Forward_Bumpato + X_T1
    S_T1_UP = F_UP + inc_T0_T1;
    S_T1_DW = F_DW + inc_T0_T1;

    % Price vanilla call with bumped forward; K1 is the vanilla strike.
    price_UP = B_T1 * mean(max(S_T1_UP - K1, 0));
    price_DW = B_T1 * mean(max(S_T1_DW - K1, 0));

elseif strcmp(exotic_name, 'PV_12M')
    % 12M forward price.
    F_base_12M = simIncrements.common.F_T0_T2;

    % Forward bump scenarios (up and down).
    F_UP = F_base_12M + bump;
    F_DW = F_base_12M - bump;

    % S_T2 = Forward_Bumpato + X_T1 + Delta_X_T1_T2
    S_T2_UP = F_UP + inc_T0_T1 + inc_T1_T2;
    S_T2_DW = F_DW + inc_T0_T1 + inc_T1_T2;

    % Price vanilla call with bumped forward; K1 is the vanilla strike.
    price_UP = B_T2 * mean(max(S_T2_UP - K1, 0));
    price_DW = B_T2 * mean(max(S_T2_DW - K1, 0));

else
    F_base = simIncrements.common.F_T0_T2;
    K2 = F_base; % Strike equals the initial ATM forward (fixed at inception).

    % Forward bump scenarios for the exotic underlying.
    F_UP = F_base + bump;
    F_DW = F_base - bump;

    % Chooser has no K1 strike.
    if strcmp(exotic_name, 'Chooser')
        K1 = [];
    end

    % Price with bumped-up forward; K2 remains fixed.
    [price_UP, ~] = price_exotic_no_grid_from_increments( ...
        exotic_name, inc_T0_T1, inc_T1_T2, ...
        B_T1, B_T2, F_UP, K2, K1, num_grid_points);

    % Price with bumped-down forward; K2 remains fixed.
    [price_DW, ~] = price_exotic_no_grid_from_increments( ...
        exotic_name, inc_T0_T1, inc_T1_T2, ...
        B_T1, B_T2, F_DW, K2, K1, num_grid_points);
end

% Central-difference delta.
delta = (price_UP - price_DW) / (2 * bump);
end