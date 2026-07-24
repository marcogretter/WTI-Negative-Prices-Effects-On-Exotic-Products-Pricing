function [alpha, beta, resultsTable] = calibrate_alpha_beta_MA(otmSurface, alpha_0, beta_0, initialGrid)
% calibrate_alpha_beta_MA
%
% Multi-start calibration of the Minimal Additive parameters.
%
% The MA calibration may be weakly identified in the common scale of
% alpha and beta. If several starting points produce essentially the
% same objective value, we select the solution closest to unit scale.

if nargin < 4 || isempty(initialGrid)
    initialGrid = [alpha_0, beta_0];
end

options = optimset('Display', 'off', 'MaxFunEvals', 1000, 'MaxIter', 1000, 'TolX', 1e-6,'TolFun', 1e-8);

nStarts = size(initialGrid, 1);

Alpha0 = initialGrid(:, 1);
Beta0  = initialGrid(:, 2);

AlphaHat = NaN(nStarts, 1);
BetaHat  = NaN(nStarts, 1);
ObjValue = NaN(nStarts, 1);
ExitFlag = NaN(nStarts, 1);

for i = 1:nStarts

    p0 = initialGrid(i, :);

    [p_results, fval, exitflag] = fminsearch(@(p) obj_fun_MA(otmSurface, p(1), p(2)), p0, options);

    AlphaHat(i) = p_results(1);
    BetaHat(i)  = p_results(2);
    ObjValue(i) = fval;
    ExitFlag(i) = exitflag;

end

AlphaBetaRatio = AlphaHat ./ BetaHat;

DistanceToUnitScale = (AlphaHat - 1).^2 + (BetaHat - 1).^2;

resultsTable = table( ...
    Alpha0, Beta0, AlphaHat, BetaHat, AlphaBetaRatio, ...
    ObjValue, ExitFlag, DistanceToUnitScale, ...
    'VariableNames', { ...
    'Alpha0', 'Beta0', ...
    'AlphaHat', 'BetaHat', 'AlphaBetaRatio', ...
    'ObjectiveValue', 'ExitFlag', 'DistanceToUnitScale' ...
    } ...
    );

minObjective = min(resultsTable.ObjectiveValue);
objectiveTolerance = 1e-8 * max(1, abs(minObjective));

isEquivalentMinimizer = abs(resultsTable.ObjectiveValue - minObjective) <= objectiveTolerance;

candidateTable = resultsTable(isEquivalentMinimizer, :);
candidateTable = sortrows(candidateTable, 'DistanceToUnitScale');

alpha = candidateTable.AlphaHat(1);
beta  = candidateTable.BetaHat(1);

resultsTable = sortrows(resultsTable, {'ObjectiveValue', 'DistanceToUnitScale'});

end