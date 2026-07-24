function point6=Risk_Management(mktData,curves,point4Params,simIncrements,vol,simParams3c,projectRoot,initialPrices)

% Point 6: Risk Management — uses the AB model for Greeks and repricing.
% Bump:
% We use absolute bumps because both the forward and the Bachelier normal
% volatility are expressed in price units, not in percentage units.
bump_vol = 0.01; % absolute bump on normal volatility, in USD/barrel/sqrt(year)
bump_fwd = 0.01; % absolute bump on forward price, in USD/barrel

front_future = curves.F_curve_all(2); % ie, first available fwd with maturity 
                                      % strictly greater than the value date 

% Choosing PV option at 6 months
% Call ATM (I choose the Call because I have more strikes of Calls than of 
% Puts on that maturity, but since buying and shorting have the same price, 
% it would have been equal to consider Puts) (or close to ATM)
Call_prices_6M = mktData.snapshot(3).C_call;
Call_strikes_6M = mktData.snapshot(3).K_call;
F_T0_T1 = curves.F_curve_all(3);

% I choose the closer Strike to the ATM one (quantitative rule)
[~, idx_ATM_6M] = min(abs(Call_strikes_6M - F_T0_T1));

% I get the strikes and the relative price of the Call
K_ATM_6M = Call_strikes_6M(idx_ATM_6M);
Prezzo_Call_ATM_6M = Call_prices_6M(idx_ATM_6M);

% PV vanilla options at 12 months
% I choose the Call for the same reason of the 12 months case
Call_prices_12M = mktData.snapshot(5).C_call;
Call_strikes_12M = mktData.snapshot(5).K_call;
F_T0_T2 = curves.F_curve_all(5);

% I choose the closer Strike to the ATM one (quantitative rule)
[~, idx_ATM_12M] = min(abs(Call_strikes_12M - F_T0_T2));

% I get the strikes and the relative price of the Call
K_ATM_12M = Call_strikes_12M(idx_ATM_12M);
Prezzo_Call_ATM_12M = Call_prices_12M(idx_ATM_12M);

% Initialization Greeks structure
Greeks = struct();

% Portfolio configuration: { 'InstrumentName', Strike, PortfolioQuantity }
% Note: K1_PoP is set via point4Params.K1_PoP.
% Negative quantity means short position.
portfolio_config = {
    'CoC',     point4Params.K1_CoC,     -1;  
    'PoP',     point4Params.K1_PoP,     -1;  
    'Chooser', [],         -1;  
    'PV_6M',   K_ATM_6M,   0;  % The PV are for hedging so initially I have
    'PV_12M',  K_ATM_12M,  0   % zero quantity of them
};

% The front future's exposure is only on the variation of the future, not
% on the volatility
Greeks.Future.delta = 1;

num_instruments = size(portfolio_config, 1);

% Pricing for all the instruments in my ptf
for i = 1:num_instruments
    instr_name = portfolio_config{i, 1};
    instr_K    = portfolio_config{i, 2};
    
    % Delta
    Greeks.(instr_name).delta = compute_delta(instr_name, simIncrements, bump_fwd, instr_K, point4Params.num_grid_points, curves.F_curve_all);
    
    % Vega 6M
    Greeks.(instr_name).vega_6M = compute_vega(instr_name, '6M', simIncrements, bump_vol, vol.AB.a_norm, vol.AB.params, simParams3c.M, simParams3c.valueSim, simParams3c.typeSim, instr_K, point4Params.num_grid_points, curves.F_curve_all);
    
    % Vega 12M (PV_6M has no Vega in 12 months)
    if strcmp(instr_name, 'PV_6M')
        Greeks.(instr_name).vega_12M = 0; 
    else
        Greeks.(instr_name).vega_12M = compute_vega(instr_name, '12M', simIncrements, bump_vol, vol.AB.a_norm, vol.AB.params, simParams3c.M, simParams3c.valueSim, simParams3c.typeSim, instr_K, point4Params.num_grid_points, curves.F_curve_all);
    end
end

% Ptf Greeks:
Greeks.Portfolio.vega_6M  = 0;
Greeks.Portfolio.vega_12M = 0;
Greeks.Portfolio.delta    = 0;

for i = 1:num_instruments
    instr_name = portfolio_config{i, 1};
    qty        = portfolio_config{i, 3};
    
    if qty ~= 0
        Greeks.Portfolio.vega_6M  = Greeks.Portfolio.vega_6M  + qty * Greeks.(instr_name).vega_6M;
        Greeks.Portfolio.vega_12M = Greeks.Portfolio.vega_12M + qty * Greeks.(instr_name).vega_12M;
        Greeks.Portfolio.delta    = Greeks.Portfolio.delta    + qty * Greeks.(instr_name).delta;
    end
end

% Cascade Hedging: 
% I start from the 12M vega and hedge it with the 12M call. Then I hedge
% the portfolio 6M vega partially with the 12M call and complete the hedge
% with the 6M call. Once the vegas are hedged, I add the portfolio delta to
% the deltas of the hedging calls and use the front future to hedge the
% remaining delta:

% Quantity to buy of the 12M call. I take the integer part because I cannot
% buy 3.1564 calls, while I can buy 3.
weight_PV_12M = - Greeks.Portfolio.vega_12M/Greeks.PV_12M.vega_12M;

if weight_PV_12M < 0
    weight_PV_12M = ceil(weight_PV_12M);
else
    weight_PV_12M = floor(weight_PV_12M);
end

% Quantity to buy of the 6 months Call to cover Vega 6M of the ptf + Vega
% 6M of the 12M PV:
weight_PV_6M = -((Greeks.Portfolio.vega_6M + weight_PV_12M*Greeks.PV_12M.vega_6M)/Greeks.PV_6M.vega_6M);

if weight_PV_6M < 0
    weight_PV_6M = ceil(weight_PV_6M);
else
    weight_PV_6M = floor(weight_PV_6M);
end

% Quantity to buy of the Front Future to cover the ptf delta + delta
% Call_12M + delta call_6M:

weight_Front_Future = - (Greeks.Portfolio.delta + ...
    weight_PV_6M*Greeks.PV_6M.delta + weight_PV_12M*Greeks.PV_12M.delta);

if weight_Front_Future < 0
    weight_Front_Future = ceil(weight_Front_Future);
else
    weight_Front_Future = floor(weight_Front_Future);
end

% Compute the hedging cost:
hedge_cost = 0.5*0.0001*abs(weight_Front_Future) + ...
    0.5*0.0004*(abs(weight_PV_12M) + abs(weight_PV_6M));

% Cash-flow on the 2nd June 2020:
CF0 = -(portfolio_config{1,3}*initialPrices.CoC_AB +...
    portfolio_config{2,3}*initialPrices.PoP_AB +...
    portfolio_config{3,3}*initialPrices.Chooser_AB +...
    weight_PV_12M*Prezzo_Call_ATM_12M + weight_PV_6M*Prezzo_Call_ATM_6M)...
    - hedge_cost;

% Testing Hedging on the next two tuesdays
% Building reset dates:
originalValuationDate = datetime(2020, 6, 2);
originalValuationDateNum = datenum(originalValuationDate);

resetDatesRaw = [
    addtodate(originalValuationDateNum, 6,  'month');
    addtodate(originalValuationDateNum, 12, 'month')
];

resetDatesNum = busdate(resetDatesRaw, 'follow');
resetDates = datetime(resetDatesNum, 'ConvertFrom', 'datenum');
firstTuesdayPrices = recalibration_new_value_date( ...
    projectRoot, datetime(2020, 6, 9),F_T0_T2, point4Params, resetDates);

% Initial prices of the exotics:
price_CoC_AB     = initialPrices.CoC_AB;
price_PoP_AB     = initialPrices.PoP_AB;
price_Chooser_AB = initialPrices.Chooser_AB;

% First tuesday: 9/6/2020
% Value of the ptf on the first tuesday:
delta_CoC_first_Tuesday = firstTuesdayPrices.price_exotic_CoC - ...
    price_CoC_AB;
delta_PoP_first_Tuesday = firstTuesdayPrices.price_exotic_PoP - ...
    price_PoP_AB;
delta_Chooser_first_Tuesday = firstTuesdayPrices.price_exotic_Chooser - ...
    price_Chooser_AB;

% delta value of the calls on the first tuesday:
[FirstTuesday_Call6M_Price, ~, ~] = get_current_call_price_for_rm(firstTuesdayPrices.mktData_new.snapshot(3), K_ATM_6M);
[FirstTuesday_Call12M_Price, ~, ~] = get_current_call_price_for_rm(firstTuesdayPrices.mktData_new.snapshot(5), K_ATM_12M);

delta_Call6M_first_Tuesday = weight_PV_6M*(FirstTuesday_Call6M_Price-Prezzo_Call_ATM_6M);
delta_Call12M_first_Tuesday = weight_PV_12M*(FirstTuesday_Call12M_Price-Prezzo_Call_ATM_12M);

% delta value of the front future:
delta_fFuture_first_Tuesday = weight_Front_Future*...
    (firstTuesdayPrices.curves_new.F_curve_all(2) - front_future);

% P&L first Tuesday, 9/6/2020
Unhedged_PnL_first = portfolio_config{1,3}*delta_CoC_first_Tuesday + ...
                     portfolio_config{2,3}*delta_PoP_first_Tuesday + ...
                     portfolio_config{3,3}*delta_Chooser_first_Tuesday;
P_L_first_Tuesday = Unhedged_PnL_first + ...
    delta_Call6M_first_Tuesday + delta_Call12M_first_Tuesday + ...
    delta_fFuture_first_Tuesday - hedge_cost;

fprintf('P&L on the 9th June 2020: %8.4f\n\n',P_L_first_Tuesday);

% P&L on the second Tuesday, ie 16th June 2020:
% Value of the ptf on the second tuesday
secondTuesdayPrices = recalibration_new_value_date( ...
    projectRoot, datetime(2020, 6, 16), F_T0_T2, point4Params, resetDates);

delta_CoC_second_Tuesday = secondTuesdayPrices.price_exotic_CoC - price_CoC_AB;
delta_PoP_second_Tuesday = secondTuesdayPrices.price_exotic_PoP - price_PoP_AB;
delta_Chooser_second_Tuesday = secondTuesdayPrices.price_exotic_Chooser - price_Chooser_AB;

% delta value of the calls on the first tuesday:
[secondTuesday_Call6M_Price, ~, ~] = get_current_call_price_for_rm(secondTuesdayPrices.mktData_new.snapshot(3), K_ATM_6M);
[secondTuesday_Call12M_Price, ~, ~] = get_current_call_price_for_rm(secondTuesdayPrices.mktData_new.snapshot(5), K_ATM_12M);

delta_Call6M_second_Tuesday = weight_PV_6M*(secondTuesday_Call6M_Price - Prezzo_Call_ATM_6M);
delta_Call12M_second_Tuesday = weight_PV_12M*(secondTuesday_Call12M_Price - Prezzo_Call_ATM_12M);

% delta value of the front future:
delta_fFuture_second_Tuesday = weight_Front_Future*...
    (secondTuesdayPrices.curves_new.F_curve_all(2) - front_future);

% P&L first Tuesday, 9/6/2020
Unhedged_PnL_second = portfolio_config{1,3}*delta_CoC_second_Tuesday + ...
                      portfolio_config{2,3}*delta_PoP_second_Tuesday + ...
                      portfolio_config{3,3}*delta_Chooser_second_Tuesday;
P_L_second_Tuesday = Unhedged_PnL_second + ...
    delta_Call6M_second_Tuesday + delta_Call12M_second_Tuesday + ...
    delta_fFuture_second_Tuesday - hedge_cost;

fprintf('P&L on the 16th June 2020: %8.4f\n\n',P_L_second_Tuesday);
%% Output structure
% Summary of hedge quantities
hedgePositions = table( ...
    ["PV_6M"; "PV_12M"; "FrontFuture"], ...
    [K_ATM_6M; K_ATM_12M; NaN], ...
    [Prezzo_Call_ATM_6M; Prezzo_Call_ATM_12M; front_future], ...
    [weight_PV_6M; weight_PV_12M; weight_Front_Future], ...
    'VariableNames', {'Instrument', 'Strike', 'InitialPrice', 'Quantity'} ...
);

% Summary of portfolio Greeks before hedging
greeksSummary = table( ...
    Greeks.Portfolio.delta, ...
    Greeks.Portfolio.vega_6M, ...
    Greeks.Portfolio.vega_12M, ...
    'VariableNames', {'PortfolioDelta', 'PortfolioVega6M', 'PortfolioVega12M'} ...
);

% Exotic price changes
exoticPnLBreakdown = table( ...
    ["CoC"; "PoP"; "Chooser"], ...
    [price_CoC_AB; price_PoP_AB; price_Chooser_AB], ...
    [firstTuesdayPrices.price_exotic_CoC; ...
     firstTuesdayPrices.price_exotic_PoP; ...
     firstTuesdayPrices.price_exotic_Chooser], ...
    [secondTuesdayPrices.price_exotic_CoC; ...
     secondTuesdayPrices.price_exotic_PoP; ...
     secondTuesdayPrices.price_exotic_Chooser], ...
    [delta_CoC_first_Tuesday; ...
     delta_PoP_first_Tuesday; ...
     delta_Chooser_first_Tuesday], ...
    [delta_CoC_second_Tuesday; ...
     delta_PoP_second_Tuesday; ...
     delta_Chooser_second_Tuesday], ...
    [portfolio_config{1,3}; portfolio_config{2,3}; portfolio_config{3,3}], ...
    'VariableNames', {'Instrument', 'InitialPrice', ...
    'Price_2020_06_09', 'Price_2020_06_16', ...
    'DeltaPrice_2020_06_09', 'DeltaPrice_2020_06_16', ...
    'PortfolioQuantity'} ...
);

% Hedging instruments price changes
hedgePnLBreakdown = table( ...
    ["PV_6M"; "PV_12M"; "FrontFuture"], ...
    [Prezzo_Call_ATM_6M; Prezzo_Call_ATM_12M; front_future], ...
    [FirstTuesday_Call6M_Price; ...
     FirstTuesday_Call12M_Price; ...
     firstTuesdayPrices.curves_new.F_curve_all(2)], ...
    [secondTuesday_Call6M_Price; ...
     secondTuesday_Call12M_Price; ...
     secondTuesdayPrices.curves_new.F_curve_all(2)], ...
    [delta_Call6M_first_Tuesday; ...
     delta_Call12M_first_Tuesday; ...
     delta_fFuture_first_Tuesday], ...
    [delta_Call6M_second_Tuesday; ...
     delta_Call12M_second_Tuesday; ...
     delta_fFuture_second_Tuesday], ...
    [weight_PV_6M; weight_PV_12M; weight_Front_Future], ...
    'VariableNames', {'Instrument', 'InitialPrice', ...
    'Price_2020_06_09', 'Price_2020_06_16', ...
    'PnL_2020_06_09', 'PnL_2020_06_16', ...
    'Quantity'} ...
);

% Final P&L summary
pnlSummary = table( ...
    ["2020-06-09"; "2020-06-16"], ...
    [Unhedged_PnL_first; Unhedged_PnL_second], ...
    [P_L_first_Tuesday; P_L_second_Tuesday], ...
    [P_L_first_Tuesday - Unhedged_PnL_first; ...
     P_L_second_Tuesday - Unhedged_PnL_second], ...
    [hedge_cost; hedge_cost], ...
    'VariableNames', {'Date', 'UnhedgedPnL', 'HedgedPnL', ...
    'HedgeImprovement', 'HedgingCost'} ...
);

% Reset-date check
resetDateCheck = table( ...
    resetDates(:), ...
    firstTuesdayPrices.tGrid(:), ...
    secondTuesdayPrices.tGrid(:), ...
    'VariableNames', {'ResetDate', 'TTM_2020_06_09', 'TTM_2020_06_16'} ...
);

% Main output
point6 = struct();

point6.Greeks = Greeks;
point6.portfolio_config = portfolio_config;

point6.weights.PV_6M = weight_PV_6M;
point6.weights.PV_12M = weight_PV_12M;
point6.weights.FrontFuture = weight_Front_Future;

point6.hedge_cost = hedge_cost;
point6.CF0 = CF0;

point6.pnl.Unhedged.firstTuesday = Unhedged_PnL_first;
point6.pnl.Unhedged.secondTuesday = Unhedged_PnL_second;
point6.pnl.Hedged.firstTuesday = P_L_first_Tuesday;
point6.pnl.Hedged.secondTuesday = P_L_second_Tuesday;

point6.priceChanges.firstTuesday.CoC = delta_CoC_first_Tuesday;
point6.priceChanges.firstTuesday.PoP = delta_PoP_first_Tuesday;
point6.priceChanges.firstTuesday.Chooser = delta_Chooser_first_Tuesday;

point6.priceChanges.secondTuesday.CoC = delta_CoC_second_Tuesday;
point6.priceChanges.secondTuesday.PoP = delta_PoP_second_Tuesday;
point6.priceChanges.secondTuesday.Chooser = delta_Chooser_second_Tuesday;

point6.tables.hedgePositions = hedgePositions;
point6.tables.greeksSummary = greeksSummary;
point6.tables.exoticPnLBreakdown = exoticPnLBreakdown;
point6.tables.hedgePnLBreakdown = hedgePnLBreakdown;
point6.tables.pnlSummary = pnlSummary;
point6.tables.resetDateCheck = resetDateCheck;

point6.repricing.firstTuesday = firstTuesdayPrices;
point6.repricing.secondTuesday = secondTuesdayPrices;
end