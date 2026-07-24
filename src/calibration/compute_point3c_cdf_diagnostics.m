function cdfDiag = compute_point3c_cdf_diagnostics(sim3c_est, sim3c_theory, simParams3c)
% Returns Lewis-FFT CDF reconstruction diagnostics for Point 3c.
% Checks probability mass outside the numerical domain and tail behaviour
% for each model (AB, GL, MA) and tail mode (Estimated, Theoretical).
% All inputs are already computed by run_point3c_simulate_increments.

    N_sim  = simParams3c.N_sim;
    models = {'AB', 'GL', 'MA'};

    cdfDiag.cdfQuality     = cdf_quality_table(sim3c_est, sim3c_theory, models);
    cdfDiag.tailComparison = tail_comparison_table(sim3c_est, sim3c_theory, models, N_sim);
end


function t = cdf_quality_table(sim3c_est, sim3c_theory, models)
% One row per (model, tail mode) pair. Estimated rows precede theoretical.

    nModels = numel(models);
    nRows   = nModels * 2;
    sims    = {sim3c_est, sim3c_theory};
    modes   = {'Estimated', 'Theoretical'};

    Model         = strings(nRows, 1);
    TailMode      = strings(nRows, 1);
    PLeft         = zeros(nRows, 1);
    PRight        = zeros(nRows, 1);
    LeftTailMass  = zeros(nRows, 1);
    RightTailMass = zeros(nRows, 1);
    TailMass      = zeros(nRows, 1);
    MinDiffCDF    = zeros(nRows, 1);
    LambdaMinus   = zeros(nRows, 1);
    LambdaPlus    = zeros(nRows, 1);
    MaxIncrement  = zeros(nRows, 1);
    Kurtosis      = zeros(nRows, 1);

    row = 0;
    for m = 1:nModels
        for k = 1:2
            row  = row + 1;
            name = models{m};
            s    = sims{k}.(name);

            Model(row)         = name;
            TailMode(row)      = modes{k};
            PLeft(row)         = s.invCDF.Pb;
            PRight(row)        = s.invCDF.Pe;
            LeftTailMass(row)  = s.invCDF.Pb;
            RightTailMass(row) = 1 - s.invCDF.Pe;
            TailMass(row)      = max(s.invCDF.Pb, 1 - s.invCDF.Pe);
            MinDiffCDF(row)    = min(diff(s.invCDF.PBlock));
            LambdaMinus(row)   = s.lambdaMinus;
            LambdaPlus(row)    = s.lambdaPlus;
            MaxIncrement(row)  = max(s.increments);
            Kurtosis(row)      = kurtosis(s.increments);
        end
    end

    t = table(Model, TailMode, PLeft, PRight, LeftTailMass, RightTailMass, TailMass, ...
              MinDiffCDF, LambdaMinus, LambdaPlus, MaxIncrement, Kurtosis, ...
              'VariableNames', {'Model','TailMode','PLeft','PRight','LeftTailMass', ...
                  'RightTailMass','TailMass','MinDiffCDF','LambdaMinus','LambdaPlus', ...
                  'MaxIncrement','Kurtosis'});
end


function t = tail_comparison_table(sim3c_est, sim3c_theory, models, N_sim)
% Compare estimated vs theoretical right-tail coefficients.
% ApproxMax = log(N_sim) / lambdaPlus gives the expected order-statistic
% maximum under pure exponential tail sampling.

    nModels              = numel(models);
    Model                = strings(nModels, 1);
    LambdaPlus_Est       = zeros(nModels, 1);
    LambdaPlus_Theory    = zeros(nModels, 1);
    TheoryOverEstimated  = zeros(nModels, 1);
    ApproxMax_EstTail    = zeros(nModels, 1);
    ApproxMax_TheoryTail = zeros(nModels, 1);

    for m = 1:nModels
        name = models{m};
        lp_e = sim3c_est.(name).lambdaPlus;
        lp_t = sim3c_theory.(name).lambdaPlus;

        Model(m)                = name;
        LambdaPlus_Est(m)       = lp_e;
        LambdaPlus_Theory(m)    = lp_t;
        TheoryOverEstimated(m)  = lp_t / lp_e;
        ApproxMax_EstTail(m)    = log(N_sim) / lp_e;
        ApproxMax_TheoryTail(m) = log(N_sim) / lp_t;
    end

    t = table(Model, LambdaPlus_Est, LambdaPlus_Theory, TheoryOverEstimated, ...
              ApproxMax_EstTail, ApproxMax_TheoryTail, ...
              'VariableNames', {'Model','LambdaPlus_Est','LambdaPlus_Theory', ...
                  'TheoryOverEstimated','ApproxMax_EstTail','ApproxMax_TheoryTail'});
end
