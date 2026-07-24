function [B, F, R2] = calibrate_single_maturity(K_call, C_mid, K_put, P_mid)
% CALIBRATE_SINGLE_MATURITY
% This function calibrates the Discount Factor (B) and Forward (F) using
% only European Call and Put prices from the derivative market. 

% Inputs:
% K_call: (vector) vector of the strikes of the Call that are available 
%         before the data pre-processing cleaning 
% C_mid:  (vector) mid Call prices related to the strikes
% K_put:  (vector) vector of the strikes of the Put that are available before the
%         data pre-processing cleaning
% P_mid:  (vector) mid Put prices related to the strikes

% Outputs:
% B:      (scalar) implicit discount factor given by the calibration
% F:      (scalar) implicit forward price given by the calibration
% R2:     (scalar) check of the R^2 according to the paper
    
% First data pre-processing liquidity criteria: 
% Survive only the options, prices and the strikes related, with price 
% greater than 0.1 "INDEX POINTS"

% Call
liq_mask_C = C_mid >= 0.1;
K_call_liq = K_call(liq_mask_C);
C_liq = C_mid(liq_mask_C);
    
% Put
liq_mask_P = P_mid >= 0.1;
K_put_liq = K_put(liq_mask_P);
P_liq = P_mid(liq_mask_P);
    
% We match the two new vectors that we got, we preserve only the prices of
% Calls and Puts that have the same strike in common
[K_clean, idx_C, idx_P] = intersect(K_call_liq, K_put_liq);
    
C_clean = C_liq(idx_C);
P_clean = P_liq(idx_P);
    
% Second data pre-processing criteria: we don't do the regression if we
% have less that 3 strikes after the liquidity criteria cleaning
if length(K_clean) < 3
    B = NaN; F = NaN; R2 = NaN;
    return;
end
    
% Compute the Synthetic Forward G(K) = C(K) - P(K)
G = C_clean - P_clean;

% Linear regression (OLS): G = -B*K + B*F
p = polyfit(K_clean, G, 1);
    
% Output
B = -p(1);
F = p(2) / B;

% We compute the R^2 in order to check that we obtain a result close to the
% one in the paper (ie, 0.99 ca)
G_fit = polyval(p, K_clean);
SSres = sum((G - G_fit).^2);
SStot = sum((G - mean(G)).^2);
R2 = 1 - (SSres / SStot);