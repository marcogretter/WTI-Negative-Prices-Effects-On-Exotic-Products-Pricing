function [x_grid, P_grid] = cdf_Lewis_FFT_from_cf(cfHandle, a, M, value, type)
% Compute CDF P(x) from characteristic function using shifted Lewis-FFT.
% Formula:
% P(x) = R_a - exp(-a*x)/(2*pi) * int exp(-i*u*x) phi(u - i*a) / (i*(u - i*a)) du

N = 2^M;

if strcmpi(type, 'x1')
    x1 = value;

    if x1 >= 0
        error('For type x1, value must be negative.');
    end

    dx = -2 * x1 / (N-1);
    du = (2*pi) / (N * dx);

elseif strcmpi(type, 'du')
    du = value;
    dx = (2*pi) / (N * du);
    x1 = -dx * (N-1) / 2;

else
    error('Invalid type. Use ''x1'', ''dx'', or ''du''.');
end

% building the x_grid (N * 1 vector) - symmetric
x1 = - dx * (N-1) / 2;
xN = - x1;
x_grid = (x1:dx:xN).';

% building the u_grid (N * 1 vector) - symmetric
u1 = - du * (N-1) / 2;
uN = - u1;
u_grid = (u1:du:uN).';

x_grid = x_grid(1:N);
u_grid = u_grid(1:N);

% compute the integrand values at u_grid - shifted by a
if a == 0
    error('The shift must be nonzero.');
end

u_shifted = u_grid - 1i * a;

G = cfHandle(u_shifted) ./ (1i .* u_shifted);
bad = ~isfinite(G) | abs(G) > 1e200;
G(bad) = 0;

% FFT approximation of the integral

n_minus_1 = (0:N-1).';

input_fft = G .* exp(-1i .* n_minus_1 * x1 * du);

Y = fft(input_fft);

integral_values = du .* exp(-1i .* u1 .* x_grid) .* Y;

% reconstruct the CDF values
if a > 0
    R_a = 1;
else
    R_a = 0;
end

P_grid = R_a - exp(-a * x_grid) ./ (2 * pi) .* real(integral_values);

P_grid = real(P_grid);

end