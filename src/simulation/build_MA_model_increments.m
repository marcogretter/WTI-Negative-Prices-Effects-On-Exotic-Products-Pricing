function modelSim = build_MA_model_increments( ...
    sigma_T1, sigma_T2, a_norm, ...
    paramsMA, tGrid_d, M, valueSim, typeSim, ...
    u_T0_T1, u_T1_T2, useTheoreticalTails)

if nargin < 11, useTheoreticalTails = false; end

%BUILD_MA_MODEL_INCREMENTS Simulate MA increments with atom management.

    modelName = 'MA';

    T0 = tGrid_d(1);
    T1 = tGrid_d(2);
    T2 = tGrid_d(3);

    sigma_T0 = 1;  % dummy value, not used when s = 0

    % CDF analyticity strip requires 0 < a < lambdaPlus. a_norm is the
    % call-pricing damping shift (negative); negating it gives +0.5*lambdaPlus.
    a_T0_T1 = -a_norm / (sigma_T1 * sqrt(T1));
    a_T1_T2 = -a_norm / (sigma_T2 * sqrt(T2));

    N_sim = numel(u_T0_T1);

    if numel(u_T1_T2) ~= N_sim
        error('u_T0_T1 and u_T1_T2 must have the same length.');
    end

    %% Bernoulli variables

    Bern_T0_T1 = ones(N_sim, 1);

    atomProb_T1_T2 = (sigma_T1^2 * T1) / (sigma_T2^2 * T2);
    jumpProb_T1_T2 = 1 - atomProb_T1_T2;

    if jumpProb_T1_T2 < -1e-12 || jumpProb_T1_T2 > 1 + 1e-12
        error('Invalid MA jump probability. Check sigma_T1, sigma_T2, T1, T2.');
    end

    jumpProb_T1_T2 = min(max(jumpProb_T1_T2, 0), 1);

    Bern_T1_T2 = double(rand(N_sim, 1) < jumpProb_T1_T2);

    %% Simulate increments

    x_T0_T1 = manage_MA_increments( ...
        Bern_T0_T1, u_T0_T1, modelName, ...
        sigma_T0, sigma_T1, a_T0_T1, ...
        paramsMA, valueSim, typeSim, M, tGrid_d(1:2), useTheoreticalTails);

    x_T1_T2 = manage_MA_increments( ...
        Bern_T1_T2, u_T1_T2, modelName, ...
        sigma_T1, sigma_T2, a_T1_T2, ...
        paramsMA, valueSim, typeSim, M, tGrid_d(2:3), useTheoreticalTails);

    %% Store

    modelSim = struct();

    modelSim.model = modelName;

    modelSim.sigma.T0 = sigma_T0;
    modelSim.sigma.T1 = sigma_T1;
    modelSim.sigma.T2 = sigma_T2;

    modelSim.a.T0_T1 = a_T0_T1;
    modelSim.a.T1_T2 = a_T1_T2;

    modelSim.bernoulli.T0_T1 = Bern_T0_T1;
    modelSim.bernoulli.T1_T2 = Bern_T1_T2;

    modelSim.prob.atom_T1_T2 = atomProb_T1_T2;
    modelSim.prob.jump_T1_T2 = jumpProb_T1_T2;

    modelSim.increments.T0_T1 = x_T0_T1(:);
    modelSim.increments.T1_T2 = x_T1_T2(:);

    modelSim.state.T1 = modelSim.increments.T0_T1;
    modelSim.state.T2 = modelSim.increments.T0_T1 + modelSim.increments.T1_T2;

end