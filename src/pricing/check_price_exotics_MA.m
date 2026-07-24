function [price_MA_numerical,price_analytic_MA] = check_price_exotics_MA(alpha_MA,beta_MA,tGrid_d,sigma_T1_MA,sigma_T2_MA,B_T1,B_T2,xtarget_MA_T0_T1,xtarget_MA_T1_T2,exoticType)
% This function checks the prices of CoC and the Chooser for some K1 in
% order 

exoticType = validatestring(exoticType, {'CoC', 'Chooser','PoP'});
% Initializing parameters:
C_MA = (alpha_MA*beta_MA)/(alpha_MA+beta_MA);
gamma_MA = (1/alpha_MA) - (1/beta_MA);
mu_T1 = gamma_MA*sigma_T1_MA*sqrt(tGrid_d(2));
mu_T2 = gamma_MA*sigma_T2_MA*sqrt(tGrid_d(3));
p_plus_T2 = beta_MA/(sigma_T2_MA*sqrt(tGrid_d(3)));
p_minus_T2 = alpha_MA/(sigma_T2_MA*sqrt(tGrid_d(3)));
p_plus_T1 = beta_MA/(sigma_T1_MA*sqrt(tGrid_d(2)));
p_minus_T1 = alpha_MA/(sigma_T1_MA*sqrt(tGrid_d(2)));
p0 = (p_minus_T2*p_plus_T2)/(p_minus_T1*p_plus_T1);
A_minus = p0*((p_minus_T1 - p_minus_T2)*(p_plus_T1 + p_minus_T2) / (p_plus_T2 + p_minus_T2));
A_plus = p0*((p_minus_T1 + p_plus_T2)*(p_plus_T1 - p_plus_T2) / (p_plus_T2 + p_minus_T2));
S = sigma_T1_MA*sqrt(tGrid_d(2));
d = mu_T2 - mu_T1;

num_grid_points = 2000; 
qMain = linspace(0.001, 0.999, round(0.80 * num_grid_points));
qTail = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01, ...
    0.99, 0.995, 0.999, 0.9995, 0.9999, 1];

switch exoticType
    case 'CoC'
        %% Numerical:
        % Evaluation of the largest K1 i can choose, price of a call
        % increases with y (ie, increment), i take the y minimum value, 
        %  y > -(mu_T2 - mu_T1) and i compute the price of the inner call 
        % in that point, at this point i select the K1 value that satisfies
        % this condition, for all the values that are below this one i am
        % ok
        
        E_1 = @(y) A_minus*(((mu_T2-mu_T1 + y)/p_minus_T2) - (1/(p_minus_T2^2)) + (1/(p_minus_T2^2))*exp(-p_minus_T2*(y+mu_T2-mu_T1))) +...
            A_plus*(((mu_T2-mu_T1+y)/p_plus_T2) + (1/(p_plus_T2^2))) + p0*(d + y);
        K1_MA = (B_T2/B_T1)*E_1(-(mu_T2-mu_T1)); 
    
        xQuantileGrid = quantile(xtarget_MA_T0_T1, unique([qMain, qTail])).';
    
        xCentralGrid = linspace( ...
            quantile(xtarget_MA_T0_T1, 0.001), ...
            quantile(xtarget_MA_T0_T1, 0.999), ...
            round(0.20 * num_grid_points)).';
    
        x_T1_grid = unique([xQuantileGrid; xCentralGrid; min(xtarget_MA_T0_T1); max(xtarget_MA_T0_T1)]);
        
        % vector to save prices in the grid
        nGrid = numel(x_T1_grid);
        
        payoff_in_T1_grid = zeros(nGrid, 1);
        
        for k = 1:nGrid
            payoff_in_T1_grid(k) = mean(max(x_T1_grid(k) + xtarget_MA_T1_T2, 0));
        end
        
        % Interpolation on the grid
        payoff_in_T1 = interp1(x_T1_grid, payoff_in_T1_grid, xtarget_MA_T0_T1, 'pchip');
        
        % Compute Call-on-Call final price
        price_MA_numerical = B_T1 * mean(max((B_T2/B_T1) .* payoff_in_T1 - K1_MA, 0));
        
        % Analytic:
        price_analytic_MA = C_MA * B_T2 * (...
            (A_plus/(p_plus_T2*beta_MA)) * exp(beta_MA*(gamma_MA + d/S)) * ((S/beta_MA) + (1/p_plus_T2)) + ...
            (A_minus/p_minus_T2) * exp(beta_MA*(gamma_MA + d/S)) * ((S/(beta_MA^2)) - (1/(p_minus_T2*beta_MA)) + 1/(p_minus_T2*(S*p_minus_T2 + beta_MA))) - ...
            (K1_MA/(beta_MA*(B_T2/B_T1))) * exp(beta_MA*(gamma_MA + d/S)) + ...
            p0 * (S/(beta_MA^2)) * exp(beta_MA*(gamma_MA + d/S)));
    case 'PoP'
        E_1 = @(y) (A_minus/(p_minus_T2^2))*exp(-p_minus_T2*(y + d));
        K1_MA = (B_T2/B_T1)*E_1(-(mu_T2-mu_T1));
        % Numerical:
        xQuantileGrid = quantile(xtarget_MA_T0_T1, unique([qMain, qTail])).';
    
        xCentralGrid = linspace( ...
            quantile(xtarget_MA_T0_T1, 0.001), ...
            quantile(xtarget_MA_T0_T1, 0.999), ...
            round(0.20 * num_grid_points)).';
    
        x_T1_grid = unique([xQuantileGrid; xCentralGrid; min(xtarget_MA_T0_T1); max(xtarget_MA_T0_T1)]);
        
        nGrid = numel(x_T1_grid);
        
        payoff_in_T1_grid = zeros(nGrid, 1);
        
        for k = 1:nGrid
            payoff_in_T1_grid(k) = mean(max(-x_T1_grid(k) - xtarget_MA_T1_T2, 0));
        end
        
        % Interpolation on the grid
        payoff_in_T1 = interp1(x_T1_grid, payoff_in_T1_grid, xtarget_MA_T0_T1, 'pchip');
        
        % Comput PoP final price
        price_MA_numerical = B_T1 * mean(max(-(B_T2/B_T1) .* payoff_in_T1 + K1_MA, 0));


        % Analytic:
        price_analytic_MA = C_MA*B_T2*exp(beta_MA*(gamma_MA + d/S))*(...
            (K1_MA/(beta_MA*(B_T2/B_T1))) - (A_minus/((p_minus_T2^2)*...
            (p_minus_T2*S + beta_MA)))*exp(beta_MA*(gamma_MA + d/S)));

    case 'Chooser'

        % Numerical 
        xQuantileGrid = quantile(xtarget_MA_T0_T1, unique([qMain, qTail])).';
    
        xCentralGrid = linspace( ...
            quantile(xtarget_MA_T0_T1, 0.001), ...
            quantile(xtarget_MA_T0_T1, 0.999), ...
            round(0.20 * num_grid_points)).';
    
        x_T1_grid = unique([xQuantileGrid; xCentralGrid; min(xtarget_MA_T0_T1); max(xtarget_MA_T0_T1)]);
        
        nGrid = numel(x_T1_grid);
        
        payoff_in_T1_grid = zeros(nGrid, 1);
        
        for k = 1:nGrid
            innerCall = mean(max((B_T1/B_T2) * x_T1_grid(k) + xtarget_MA_T1_T2, 0));
            innerPut  = mean(max(-(B_T1/B_T2) * x_T1_grid(k) - xtarget_MA_T1_T2, 0));
        
            payoff_in_T1_grid(k) = max(innerCall, innerPut);
        end
        % Interpolation on the grid
        payoff_in_T1_chooser = interp1(x_T1_grid, payoff_in_T1_grid, xtarget_MA_T0_T1, 'pchip');
        
        price_MA_numerical = B_T2 * mean(payoff_in_T1_chooser);


        % Analytic
        price_analytic_MA = ...
        (C_MA/(beta_MA^2))*exp(beta_MA*gamma_MA)*(B_T2*sigma_T2_MA*sqrt(tGrid_d(3))+...
        B_T1*sigma_T1_MA*sqrt(tGrid_d(2)));
end
