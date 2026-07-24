function calibDiag = compute_point2_calibration_diagnostics(vol, fftParams)
% Returns calibration diagnostics for AB, MA and GL after point 2.
% Does not modify calibration results or model parameters.
% Residuals: model call-equivalent price minus market call-equivalent price.

[rAB, rMA, rGL] = model_residuals(vol, fftParams);
otm = vol.otmSurface;

calibDiag.global      = global_error_table(rAB, rMA, rGL);
calibDiag.byMaturity  = group_error_table('MaturityCode', maturity_labels(otm),       rAB, rMA, rGL);
calibDiag.byChiBucket = group_error_table('ChiBucket',    chi_bucket_labels(otm.chi), rAB, rMA, rGL);
calibDiag.bySide      = group_error_table('Side',         side_labels(otm.IsCall),    rAB, rMA, rGL);
if isfield(vol.MA, 'multistartTable') && ~isempty(vol.MA.multistartTable)
    calibDiag.maMultistartTable = multistart_report_table(vol.MA.multistartTable);
end
end

function [rAB, rMA, rGL] = model_residuals(vol, fftParams)
% Recompute model prices using the same formula as the calibration objectives,
% then return residuals (model price minus market call-equivalent price).
% Each model price: B * sigmaATM * sqrt(T) * c(I0 * chi) / I0.

otm   = vol.otmSurface;
scale = otm.DiscountFactor .* otm.sigmaATM .* sqrt(otm.TTM);
M = fftParams.M;  value = fftParams.value;  type = fftParams.type;

c_AB = price_Bachelier_Lewis_FFT(vol.AB.I0 .* otm.chi, 1, 1, vol.AB.a_norm, 1, ...
    vol.AB.eta, vol.AB.kappa, M, value, type);
pAB  = scale .* c_AB ./ vol.AB.I0;

c_MA = call_MA_normalized(vol.MA.I0 .* otm.chi, vol.MA.alpha, vol.MA.beta);
pMA  = scale .* c_MA ./ vol.MA.I0;

c_GL = price_Bachelier_Lewis_FFT(vol.GL.I0 .* otm.chi, 1, 1, vol.GL.a_norm, 1, ...
    vol.GL.alpha, vol.GL.beta, M, value, type, 'GL');
pGL  = scale .* c_GL ./ vol.GL.I0;

rAB = pAB - otm.allCall;
rMA = pMA - otm.allCall;
rGL = pGL - otm.allCall;
end

function t = global_error_table(rAB, rMA, rGL)
% Global fit quality: one row per model across the entire OTM surface.

N = numel(rAB);
t = table( ...
    {'AB'; 'MA'; 'GL'}, repmat(N, 3, 1), ...
    [compute_rmse(rAB); compute_rmse(rMA); compute_rmse(rGL)], ...
    [compute_mae(rAB);  compute_mae(rMA);  compute_mae(rGL)], ...
    [max(abs(rAB)); max(abs(rMA)); max(abs(rGL))], ...
    [mean(rAB); mean(rMA); mean(rGL)], ...
    [median(rAB); median(rMA); median(rGL)], ...
    [mean(rAB > 0); mean(rMA > 0); mean(rGL > 0)], ...
    'VariableNames', {'Model','N','RMSE_Price','MAE_Price', ...
    'MaxAbsError','MeanError','MedianError','PositiveErrorShare'});
end

function t = group_error_table(groupVarName, groupLabels, rAB, rMA, rGL)
% Fit quality broken down by any grouping variable.
% Rows appear in order of first occurrence of each group label.

groups = unique(groupLabels, 'stable');
nG     = numel(groups);

GroupLabel   = groups;
N            = zeros(nG, 1);
RMSE_AB      = zeros(nG, 1);  RMSE_MA      = zeros(nG, 1);  RMSE_GL      = zeros(nG, 1);
MAE_AB       = zeros(nG, 1);  MAE_MA       = zeros(nG, 1);  MAE_GL       = zeros(nG, 1);
MeanError_AB = zeros(nG, 1);  MeanError_MA = zeros(nG, 1);  MeanError_GL = zeros(nG, 1);

for i = 1:nG
    mask = groupLabels == groups(i);
    N(i)            = sum(mask);
    RMSE_AB(i)      = compute_rmse(rAB(mask));  RMSE_MA(i)      = compute_rmse(rMA(mask));  RMSE_GL(i)      = compute_rmse(rGL(mask));
    MAE_AB(i)       = compute_mae(rAB(mask));   MAE_MA(i)       = compute_mae(rMA(mask));   MAE_GL(i)       = compute_mae(rGL(mask));
    MeanError_AB(i) = mean(rAB(mask));          MeanError_MA(i) = mean(rMA(mask));          MeanError_GL(i) = mean(rGL(mask));
end

t = table(GroupLabel, N, RMSE_AB, RMSE_MA, RMSE_GL, MAE_AB, MAE_MA, MAE_GL, ...
    MeanError_AB, MeanError_MA, MeanError_GL, ...
    'VariableNames', {groupVarName, 'N', ...
    'RMSE_AB','RMSE_MA','RMSE_GL', ...
    'MAE_AB','MAE_MA','MAE_GL', ...
    'MeanError_AB','MeanError_MA','MeanError_GL'});
end

function labels = maturity_labels(otm)
% Map each surface row to the expiry date string of its maturity.

mats   = unique(otm.MaturityIndex, 'stable');
labels = strings(numel(otm.MaturityIndex), 1);
for j = 1:numel(mats)
    mask         = otm.MaturityIndex == mats(j);
    labels(mask) = string(datestr(otm.Expiry(find(mask, 1)), 'dd-mmm-yyyy'));
end
end


function labels = chi_bucket_labels(chi)
% Assign each normalised moneyness value to a moneyness bucket.

edges  = [-Inf, -2, -1, 0, 1, 2, Inf];
names  = {'chi<-2', '-2<=chi<-1', '-1<=chi<0', '0<=chi<1', '1<=chi<2', 'chi>=2'};
labels = strings(numel(chi), 1);
for b  = 1:numel(names)
    labels(chi >= edges(b) & chi < edges(b+1)) = names{b};
end
end


function labels = side_labels(isCall)
% Convert the IsCall flag to a readable string per row.

isCall          = logical(isCall);
labels          = strings(numel(isCall), 1);
labels(isCall)  = 'Call';
labels(~isCall) = 'Put';
end

function t = multistart_report_table(ms)
% Return one row per MA starting point with calibrated parameters and objective.

ms = sortrows(ms, {'Alpha0', 'Beta0'});

t = table( ...
    ms.Alpha0, ...
    ms.Beta0, ...
    ms.AlphaHat, ...
    ms.BetaHat, ...
    ms.AlphaBetaRatio, ...
    ms.ObjectiveValue, ...
    'VariableNames', { ...
    'Alpha0', ...
    'Beta0', ...
    'AlphaStar', ...
    'BetaStar', ...
    'AlphaStarOverBetaStar', ...
    'ObjectiveValue' ...
    } ...
    );
end

function v = compute_rmse(r)
    v = sqrt(mean(r .^ 2));
end

function v = compute_mae(r)
    v = mean(abs(r));
end