function modelSim = build_fft_model_increments(modelName, sigma_T1, sigma_T2, a_norm, params, tGrid_d, M, valueSim, typeSim, u_T0_T1, u_T1_T2, useTheoreticalTails)

if nargin < 12, useTheoreticalTails = false; end

%BUILD_FFT_MODEL_INCREMENTS Simulate additive increments for AB/GL models.
%
% Output fields:
%   modelSim.model
%   modelSim.sigma.T0
%   modelSim.sigma.T1
%   modelSim.sigma.T2
%   modelSim.a.T0_T1
%   modelSim.a.T1_T2
%   modelSim.grid.T0_T1.x
%   modelSim.grid.T0_T1.P
%   modelSim.grid.T1_T2.x
%   modelSim.grid.T1_T2.P
%   modelSim.increments.T0_T1
%   modelSim.increments.T1_T2

    modelName = upper(string(modelName));

    sigma_T0 = 1;  % dummy value, not used when s = 0

    T0 = tGrid_d(1);
    T1 = tGrid_d(2);
    T2 = tGrid_d(3);

    % CDF analyticity strip requires 0 < a < lambdaPlus. a_norm is the
    % call-pricing damping shift (negative); negating it gives +0.5*lambdaPlus.
    a_T0_T1 = -a_norm / (sigma_T1 * sqrt(T1));
    a_T1_T2 = -a_norm / (sigma_T2 * sqrt(T2));

    % Increment T0 -> T1

    cfForFFT_T0_T1 = @(u) cf_increment_model(modelName, u, T0, T1, sigma_T0, sigma_T1, params).phiFFT;

    [x_grid_T0_T1, P_grid_T0_T1] = cdf_Lewis_FFT_from_cf(cfForFFT_T0_T1, a_T0_T1, M, valueSim, typeSim);

    if useTheoreticalTails
        [lm, lp] = fft_theoretical_lambdas(modelName, params, sigma_T1, T1);
        invCDF_T0_T1 = build_inverse_cdf_spline(x_grid_T0_T1, P_grid_T0_T1, lm, lp);
    else
        invCDF_T0_T1 = build_inverse_cdf_spline(x_grid_T0_T1, P_grid_T0_T1);
    end

    x_T0_T1 = simulated_increments(u_T0_T1, invCDF_T0_T1);

    % Increment T1 -> T2

    cfForFFT_T1_T2 = @(u) cf_increment_model(modelName, u, T1, T2, sigma_T1, sigma_T2, params).phiFFT;

    [x_grid_T1_T2, P_grid_T1_T2] = cdf_Lewis_FFT_from_cf(cfForFFT_T1_T2, a_T1_T2, M, valueSim, typeSim);

    if useTheoreticalTails
        [lm, lp] = fft_theoretical_lambdas(modelName, params, sigma_T2, T2);
        invCDF_T1_T2 = build_inverse_cdf_spline(x_grid_T1_T2, P_grid_T1_T2, lm, lp);
    else
        invCDF_T1_T2 = build_inverse_cdf_spline(x_grid_T1_T2, P_grid_T1_T2);
    end

    x_T1_T2 = simulated_increments(u_T1_T2, invCDF_T1_T2);

    % Store

    modelSim = struct();

    modelSim.model = char(modelName);

    modelSim.sigma.T0 = sigma_T0;
    modelSim.sigma.T1 = sigma_T1;
    modelSim.sigma.T2 = sigma_T2;

    modelSim.a.T0_T1 = a_T0_T1;
    modelSim.a.T1_T2 = a_T1_T2;

    modelSim.grid.T0_T1.x = x_grid_T0_T1;
    modelSim.grid.T0_T1.P = P_grid_T0_T1;

    modelSim.grid.T1_T2.x = x_grid_T1_T2;
    modelSim.grid.T1_T2.P = P_grid_T1_T2;

    modelSim.invCDF.T0_T1 = invCDF_T0_T1;
    modelSim.invCDF.T1_T2 = invCDF_T1_T2;

    modelSim.increments.T0_T1 = x_T0_T1(:);
    modelSim.increments.T1_T2 = x_T1_T2(:);

    modelSim.state.T1 = modelSim.increments.T0_T1;
    modelSim.state.T2 = modelSim.increments.T0_T1 + modelSim.increments.T1_T2;

end

function [lm, lp] = fft_theoretical_lambdas(modelName, params, sigma, t)
    qt = sigma * sqrt(t);
    if strcmpi(modelName, 'AB')
        sq = sqrt(params.eta^2 + 1/params.kappa);
        lm = (-params.eta + sq) / qt;
        lp = ( params.eta + sq) / qt;
    else 
        lm = params.alpha / qt;
        lp = params.beta  / qt;
    end
end