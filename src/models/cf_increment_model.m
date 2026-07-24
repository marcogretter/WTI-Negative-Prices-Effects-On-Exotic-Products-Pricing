function inc = cf_increment_model(model, u, s, t, sigma_s, sigma_t, params)
% Build increment CF phi_{s,t}(u) for AB, MA, GL.
% inc.phiTotal : total increment CF
% inc.phiFFT   : CF to pass to FFT
% inc.hasAtom  : true if atom must be simulated separately
% inc.p0       : atom probability

% In the MA case we remove the atom before FFT.

model = upper(model);

switch model

    case 'AB'
        eta   = params.eta;
        kappa = params.kappa;

        logphi_t = additiveBachelierLogCf(u, t, sigma_t, eta, kappa);

        if s == 0
            logphi_s = zeros(size(u));
        else
            logphi_s = additiveBachelierLogCf(u, s, sigma_s, eta, kappa);
        end

        phiInc = exp(logphi_t - logphi_s);

        inc.phiTotal = phiInc;
        inc.phiFFT   = phiInc;
        inc.hasAtom  = false;
        inc.p0       = 0;


    case 'MA'
        alpha = params.alpha;
        beta  = params.beta;

        phi_t = cf_MA_marginal(u, t, sigma_t, alpha, beta);

        if s == 0
            phi_s = ones(size(u));
        else
            phi_s = cf_MA_marginal(u, s, sigma_s, alpha, beta);
        end

        phiInc = phi_t ./ phi_s;

        q_s = sigma_s * sqrt(s); % scale for s
        q_t = sigma_t * sqrt(t); % scale for t

        if s == 0 || q_s == 0
            p0 = 0;
        else
            p0 = (q_s / q_t)^2;
        end

        gammaLoc = 1/alpha - 1/beta;

        x0 = gammaLoc * (q_t - q_s); % location of the atom

        if p0 > 0
            atomCF = p0 .* exp(1i .* u .* x0);

            phiFFT = (phiInc - atomCF) ./ (1 - p0); % absolutely continuous part for FFT

            hasAtom = true;
        else
            phiFFT = phiInc;
            x0 = 0;
            hasAtom = false;
        end

        inc.phiTotal = phiInc;
        inc.phiFFT   = phiFFT;
        inc.hasAtom  = hasAtom;
        inc.p0       = p0;
        inc.x0       = x0;


    case 'GL'
        alpha = params.alpha;
        beta  = params.beta;

        logphi_t = cf_GL_marginal_log(u, t, sigma_t, alpha, beta);

        if s == 0
            logphi_s = zeros(size(u));
        else
            logphi_s = cf_GL_marginal_log(u, s, sigma_s, alpha, beta);
        end

        phiInc = exp(logphi_t - logphi_s);

        inc.phiTotal = phiInc;
        inc.phiFFT   = phiInc;
        inc.hasAtom  = false;
        inc.p0       = 0;


    otherwise
        error('Unknown model. Use AB, MA, or GL.');
end

end
