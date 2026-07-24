function [price_exotic, diagnostics] = price_exotic_mesh_from_increments( ...
    exoticType, x_T0_T1, x_T1_T2, ...
    B_T1, B_T2, F_T0_T2, K2, K1, num_grid_points)

%PRICE_EXOTIC_MESH_FROM_INCREMENTS Numerical stochastic-mesh pricing
% for CoC, PoP and Chooser from already simulated additive increments.
%
% Assumption:
%   The simulated process is linear:
%       X_T2 = X_T1 + DeltaX_T1_T2
%
%   The terminal underlying is:
%       S_T2 = F(T0,T2) + X_T2
%
% Therefore:
%       S_T2 - K2 = X_T2 + F(T0,T2) - K2
%
% Inputs:
%   exoticType      : 'CoC', 'PoP', or 'Chooser'
%   x_T0_T1         : N_sim x 1 simulated increment from T0 to T1
%   x_T1_T2         : N_sim x 1 simulated increment from T1 to T2
%   B_T1, B_T2      : discount factors B(0,T1), B(0,T2)
%   F_T0_T2         : forward F(T0,T2)
%   K2              : strike of the inner vanilla option
%   K1              : strike of the compound option, ignored for Chooser
%   num_grid_points : number of grid points for stochastic mesh
%
% Outputs:
%   price_exotic    : numerical exotic price
%   diagnostics     : struct with useful checks

    if nargin < 9 || isempty(num_grid_points)
        num_grid_points = 2000;
    end

    if nargin < 8 || isempty(K1)
        K1 = 0;
    end

    exoticType = validatestring(exoticType, {'CoC', 'PoP', 'Chooser'});

    x_T1 = x_T0_T1(:);
    dx_T1_T2 = x_T1_T2(:);

    if numel(x_T1) ~= numel(dx_T1_T2)
        error('x_T0_T1 and x_T1_T2 must have the same length.');
    end

    if B_T1 <= 0 || B_T2 <= 0
        error('Discount factors must be positive.');
    end

    if num_grid_points < 10
        error('num_grid_points is too small.');
    end

    B12 = B_T2 / B_T1;

    % Moneyness adjustment:
    % S_T2 - K2 = X_T2 + F_T0_T2 - K2
    m2 = F_T0_T2 - K2;

    % Build stochastic mesh on X_T1
    
    qMain = linspace(0.001, 0.999, round(0.80 * num_grid_points));
    qTail = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, ...
             0.99, 0.995, 0.999, 0.9995, 0.9999, 1];
    
    xQuantileGrid = quantile(x_T1, unique([qMain, qTail])).';
    
    xCentralGrid = linspace( ...
        quantile(x_T1, 0.001), ...
        quantile(x_T1, 0.999), ...
        round(0.20 * num_grid_points)).';
    
    x_T1_grid = unique([xQuantileGrid; xCentralGrid; min(x_T1); max(x_T1)]);
    num_grid_points_effective = numel(x_T1_grid);
    
    call_in_T1_grid = zeros(num_grid_points_effective, 1);
    put_in_T1_grid  = zeros(num_grid_points_effective, 1);
    
    for k = 1:num_grid_points_effective
        X_T2_cond = x_T1_grid(k) + dx_T1_T2;
    
        call_in_T1_grid(k) = mean(max( X_T2_cond + m2, 0));
        put_in_T1_grid(k)  = mean(max(-X_T2_cond - m2, 0));
    end
    
    call_in_T1 = interp1(x_T1_grid, call_in_T1_grid, x_T1, 'pchip');
    put_in_T1  = interp1(x_T1_grid, put_in_T1_grid,  x_T1, 'pchip');


% Price selected exotic and MC standard error
%
% We explicitly build the Monte Carlo variable Y_i.
% Then:
%   price = discountExotic * mean(Y_i)
%   SE    = discountExotic * std(Y_i) / sqrt(N)

    switch exoticType
    
        case 'CoC'
            % Y_i^CoC = [B(T1,T2) * C_T1^{(i)}(K2) - K1]^+
            Y_exotic = max(B12 .* call_in_T1 - K1, 0);
            discountExotic = B_T1;
    
        case 'PoP'
            % Y_i^PoP = [K1 - B(T1,T2) * P_T1^{(i)}(K2)]^+
            Y_exotic = max(K1 - B12 .* put_in_T1, 0);
            discountExotic = B_T1;
    
        case 'Chooser'
            % Y_i^Chooser = max(C_T1^{(i)}(K2), P_T1^{(i)}(K2))
            Y_exotic = max(call_in_T1, put_in_T1);
            discountExotic = B_T2;
    end
    
    N_outer = numel(Y_exotic);
    
    price_exotic = discountExotic * mean(Y_exotic);
    
    mcStdError = discountExotic * std(Y_exotic, 0) / sqrt(N_outer);
    
    ci95HalfWidth = 1.96 * mcStdError;
    ci95Low  = price_exotic - ci95HalfWidth;
    ci95High = price_exotic + ci95HalfWidth;
    % Diagnostics

    X_T2_from_steps = x_T1 + dx_T1_T2;

    diagnostics = struct();

    diagnostics.exoticType = exoticType;
    diagnostics.K1 = K1;
    diagnostics.K2 = K2;
    diagnostics.F_T0_T2 = F_T0_T2;
    diagnostics.m2 = m2;

    diagnostics.N_sim = numel(x_T1);
    diagnostics.num_grid_points = num_grid_points;
    diagnostics.num_grid_points_effective = num_grid_points_effective;

    diagnostics.Y_exotic_mean = mean(Y_exotic);
    diagnostics.Y_exotic_std  = std(Y_exotic, 0);

    diagnostics.discountExotic = discountExotic;
    
    diagnostics.mcStdError = mcStdError;
    diagnostics.ci95_half_width = ci95HalfWidth;
    diagnostics.ci95_low = ci95Low;
    diagnostics.ci95_high = ci95High;
    
    diagnostics.priceMean = price_exotic;

    diagnostics.call_in_T1_min = min(call_in_T1);
    diagnostics.call_in_T1_max = max(call_in_T1);
    diagnostics.put_in_T1_min = min(put_in_T1);
    diagnostics.put_in_T1_max = max(put_in_T1);

    diagnostics.vanilla_call_T2_direct_MC = B_T2 * mean(max( X_T2_from_steps + m2, 0));
    diagnostics.vanilla_put_T2_direct_MC  = B_T2 * mean(max(-X_T2_from_steps - m2, 0));

    % Useful identity when K2 = F_T0_T2
    diagnostics.chooser_identity_MC = B_T2 * ( ...
        mean(max(X_T2_from_steps + m2, 0)) + ...
        mean(max(-x_T1 - m2, 0)) );

end