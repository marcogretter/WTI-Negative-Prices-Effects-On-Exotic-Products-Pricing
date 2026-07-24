function [B_curve, F_curve, R2_curve] = calibrate_discount_curve(K_call_cells, C_mid_cells, K_put_cells, P_mid_cells)
% CALIBRATE_DISCOUNT_CURVE
% This function iterates on different maturities in order to calibrate the entire curve 
% of Discount Factor and Forwards.
%
% Inputs:
% K_call_cells:     (cell array) in every cell there is the vector of
%                   strikes of the Calls related to a particular maturity.
% C_mid_cells:      (cell array) in every cell there is the vector of
%                   prices of the Calls related to a particular maturity.
% K_put_cells:      (cell array) in every cell there is the vector of
%                   strikes of the Puts related to a particular maturity.
% P_mid_cells:      (cell array) in every cell there is the vector of
%                   prices of the Puts related to a particular maturity.
% Outputs:
% B_curve:          (vector) vector of the implied discounts related to the
%                   maturity
% F_curve:          (vector) vector of the implied forwards related to the
%                   maturity
% R2_curve:         (vector) vector of the R^2 for each regression, this
%                   output is used as a check

% Total number of maturities that we have:
num_maturities = length(K_call_cells);
    
% Initialize the output vectors before the for cycle
B_curve  = NaN(num_maturities, 1);
F_curve  = NaN(num_maturities, 1);
R2_curve = NaN(num_maturities, 1);
    
% We regress the implied discounts and forward at each maturity
for i = 1:num_maturities
    % Vectors of strikes and prices for the maturity
    K_c = K_call_cells{i};
    C_m = C_mid_cells{i};
    K_p = K_put_cells{i};
    P_m = P_mid_cells{i};
        
    % Calibrate B, F and check the R^2
    [B, F, R2] = calibrate_single_maturity(K_c, C_m, K_p, P_m);
    
    % Save the results:
    B_curve(i)  = B;
    F_curve(i)  = F;
    R2_curve(i) = R2;
end