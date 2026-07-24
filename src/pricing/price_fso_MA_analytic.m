function V0 = price_fso_MA_analytic(F0, K2, B02, B12, mu1, mu2, p1m, p1p, p2m, p2p)
%PRICE_FSO_MA_ANALYTIC Analytical price of a forward-start option
% in the Minimal Additive model.
%
% Payoff at T2:
%
%   [ S_T2 - K2 F(T1,T2) ]^+
%
% Therefore:
%
%   payoff = [ X + (1-K2)F0 + (1-K2/B12)Y ]^+
%
% where:
%
%   Y = f_T1
%   X = f_T2 - f_T1
%
% and X is independent of Y.
%
% Inputs:
%   F0   = F(0,T2)
%   K2   = forward-start strike multiplier
%   B02  = discount factor B(0,T2)
%   B12  = discount factor B(T1,T2)
%   mu1  = location of f_T1
%   mu2  = location of f_T2
%   p1m  = p_1^- tail parameter of f_T1
%   p1p  = p_1^+ tail parameter of f_T1
%   p2m  = p_2^- tail parameter of f_T2 / increment CDF
%   p2p  = p_2^+ tail parameter of f_T2 / increment CDF
%
% Output:
%   V0   = analytical time-0 price
%
% Increment CDF:
%
%   F_X(x) = A exp(p2m (x-m)),             x < m
%          = 1 - B,                        x = m
%          = 1 - B exp(-p2p (x-m)),        x > m
%
%   where m = mu2 - mu1.
%
% The atom at m is 1-A-B.

    tol = 1e-12;

    % Basic quantities

    h = K2 - 1.0;
    g = K2 / B12 - 1.0;

    m = mu2 - mu1;

    % Coefficients A and B of the increment CDF

    a = p2m / p1m;
    b = p2p / p1p;

    q = p2p / (p2m + p2p);
    r = p2m / (p2m + p2p);

    A = (1.0 - a) * ( b + (1.0 - b) * q );

    B = (1.0 - b) * ( a + (1.0 - a) * r );

    atom = 1.0 - A - B;

    if atom < -1e-10
        warning('Atom mass 1-A-B is negative. Check monotonicity of tail parameters.');
    end

    % Call on increment C_X(k)

    CX = @(k) call_increment(k, m, p2m, p2p, A, B);

    % Special case g = 0, i.e. K2 = B12

    if abs(g) < tol
        k0 = h * F0;
        V0 = B02 * CX(k0);
        return;
    end

    % Distribution of Y = f_T1

    C1 = 1.0 / (1.0 / p1m + 1.0 / p1p);

    DL = C1 * exp(-p1m * mu1);
    DR = C1 * exp( p1p * mu1);

    % Switching point y*

    yStar = (m - h * F0) / g;

    % Coefficients of C_X(k(y))

    M0 = m - h * F0 - A / p2m + B / p2p;

    M2 = (A / p2m) * exp(p2m * (h * F0 - m));

    N2 = (B / p2p) * exp(p2p * m - p2p * h * F0);

    % Primitive functions

    PhiLL = @(y) primitive_LL(y, DL, M0, M2, g, p1m, p2m, tol);
    PhiLR = @(y) primitive_LR(y, DR, M0, M2, g, p1p, p2m, tol);
    PhiRL = @(y) primitive_RL(y, DL, N2, g, p1m, p2p, tol);
    PhiRR = @(y) primitive_RR(y, DR, N2, g, p1p, p2p, tol);

    % Case distinction

    if g > 0

        if yStar <= mu1
            % Case 1:
            % (-inf,y*)     : call-left,  density-left
            % (y*,mu1)      : call-right, density-left
            % (mu1,+inf)    : call-right, density-right

            Pi = PhiLL(yStar) + PhiRL(mu1) - PhiRL(yStar) - PhiRR(mu1);

        else
            % Case 2:
            % (-inf,mu1)    : call-left,  density-left
            % (mu1,y*)      : call-left,  density-right
            % (y*,+inf)     : call-right, density-right

            Pi = PhiLL(mu1) + PhiLR(yStar) - PhiLR(mu1) - PhiRR(yStar);
        end

    else

        if yStar <= mu1
            % Case 3:
            % (-inf,y*)     : call-right, density-left
            % (y*,mu1)      : call-left,  density-left
            % (mu1,+inf)    : call-left,  density-right

            Pi = PhiRL(yStar) + PhiLL(mu1) - PhiLL(yStar) - PhiLR(mu1);

        else
            % Case 4:
            % (-inf,mu1)    : call-right, density-left
            % (mu1,y*)      : call-right, density-right
            % (y*,+inf)     : call-left,  density-right

            Pi = PhiRL(mu1) + PhiRR(yStar) - PhiRR(mu1) - PhiLR(yStar);
        end
    end

    V0 = B02 * Pi;

end


% Local function: call on increment

function C = call_increment(k, m, pm, pp, A, B)
%CALL_INCREMENT Computes C_X(k) = E[(X-k)^+]

    if k < m
        C = (m - k) - (A / pm) * (1.0 - exp(pm * (k - m))) + B / pp;
    else
        C = (B / pp) * exp(-pp * (k - m));
    end

end


% Local primitive Phi_LL
% call-left, density-left

function val = primitive_LL(y, DL, M0, M2, g, p1m, p2m, tol)

    lambda = p1m + p2m * g;

    term1 = M0 * exp(p1m * y) / p1m;

    term2 = -g * exp(p1m * y) * ( y / p1m - 1.0 / (p1m^2) );

    if abs(lambda) < tol
        term3 = M2 * y;
    else
        term3 = M2 * exp(lambda * y) / lambda;
    end

    val = DL * (term1 + term2 + term3);

end


% Local primitive Phi_LR
% call-left, density-right

function val = primitive_LR(y, DR, M0, M2, g, p1p, p2m, tol)

    lambda = p2m * g - p1p;

    term1 = -M0 * exp(-p1p * y) / p1p;

    term2 = g * exp(-p1p * y) * ( y / p1p + 1.0 / (p1p^2) );

    if abs(lambda) < tol
        term3 = M2 * y;
    else
        term3 = M2 * exp(lambda * y) / lambda;
    end

    val = DR * (term1 + term2 + term3);

end


% Local primitive Phi_RL
% call-right, density-left

function val = primitive_RL(y, DL, N2, g, p1m, p2p, tol)

    lambda = p1m - p2p * g;

    if abs(lambda) < tol
        val = DL * N2 * y;
    else
        val = DL * N2 * exp(lambda * y) / lambda;
    end

end


% Local primitive Phi_RR
% call-right, density-right

function val = primitive_RR(y, DR, N2, g, p1p, p2p, tol)

    lambda = p1p + p2p * g;

    if abs(lambda) < tol
        val = DR * N2 * y;
    else
        val = -DR * N2 * exp(-lambda * y) / lambda;
    end

end