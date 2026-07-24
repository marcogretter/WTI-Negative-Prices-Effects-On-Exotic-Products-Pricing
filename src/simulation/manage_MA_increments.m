function [x_target_MA] = manage_MA_increments(B, u, model_MA, ...
    sigma_s_MA, sigma_t_MA, a_MA_sim, paramsMA, ...
    valueSim, typeSim, M, tGrid, useTheoreticalTails)
%MANAGE_MA_INCREMENTS Simulates MA increments.
%
% If useTheoreticalTails = true, the inverse CDF spline uses the
% theoretical MA asymptotic tail coefficients:
%
%   lambdaMinus = alpha / (sigma_t * sqrt(t))
%   lambdaPlus  = beta  / (sigma_t * sqrt(t))
%
% Otherwise, build_inverse_cdf_spline estimates the tail coefficients
% from the FFT CDF grid.

    if nargin < 12
        useTheoreticalTails = false;
    end

    cfForFFT = @(u) cf_increment_model( ...
        model_MA, u, ...
        tGrid(1), tGrid(2), ...
        sigma_s_MA, sigma_t_MA, ...
        paramsMA ...
    ).phiFFT;

    [x_grid_MA, P_grid_MA] = cdf_Lewis_FFT_from_cf( ...
        cfForFFT, a_MA_sim, M, valueSim, typeSim ...
    );

    if useTheoreticalTails
        q_t = sigma_t_MA * sqrt(tGrid(2));

        lambdaMinus_MA = paramsMA.alpha / q_t;
        lambdaPlus_MA  = paramsMA.beta  / q_t;

        invCDF_MA = build_inverse_cdf_spline( ...
            x_grid_MA, P_grid_MA, ...
            lambdaMinus_MA, lambdaPlus_MA ...
        );
    else
        invCDF_MA = build_inverse_cdf_spline(x_grid_MA, P_grid_MA);
    end

    if tGrid(1) == 0

        x_target_MA = simulated_increments(u, invCDF_MA);

    else

        atomValue = (1 / paramsMA.alpha - 1 / paramsMA.beta) * ...
            (sigma_t_MA * sqrt(tGrid(2)) - sigma_s_MA * sqrt(tGrid(1)));

        % Start from the atom value everywhere
        x_target_MA = atomValue * ones(size(u));

        % Only where B == 1, simulate from the absolutely continuous part
        idxContinuous = (B == 1);

        if any(idxContinuous)
            x_target_MA(idxContinuous) = simulated_increments( ...
                u(idxContinuous), invCDF_MA ...
            );
        end

    end

end