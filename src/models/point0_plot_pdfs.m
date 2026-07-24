function point0 = point0_plot_pdfs()
% Point 0: asymptotic tail coefficients and log-linear PDF comparison.
% Evaluates MA, GL, and Additive Bachelier PDFs with common p+ = 1.5, p- = 0.9.

pPlus  = 1.5;
pMinus = 0.9;

scale = 1.0;

alphaTail = pMinus * scale;
betaTail  = pPlus  * scale;

% Additive Bachelier parameters, alpha = 1/2
etaAB = (pPlus - pMinus) / 2;
kAB   = 1 / (pPlus * pMinus);

x = linspace(-12, 12, 1000).';

% MA and GL explicit pdfs
fMA = pdfMA(x, alphaTail, betaTail);
fGL = pdfGL(x, alphaTail, betaTail);

% AB pdf via Fourier inversion / FFT
M = 15;
value = -20;      % left endpoint of x/y grid
type = 'x1';

fAB = pdf_Bachelier_FFT(x, etaAB, kAB, M, value, type);

fig = figure;
semilogy(x, fMA, 'LineWidth', 1.4); hold on;
semilogy(x, fGL, 'LineWidth', 1.4);
semilogy(x, fAB, 'LineWidth', 1.4);
grid on;

xlabel('x');
ylabel('pdf');
title('Log-linear PDF comparison');
legend('MA', 'GL', 'Additive Bachelier', 'Location', 'best');

point0.tailCoefficients.pPlus  = pPlus;
point0.tailCoefficients.pMinus = pMinus;

point0.parameters.pPlus     = pPlus;
point0.parameters.pMinus    = pMinus;
point0.parameters.scale     = scale;
point0.parameters.alphaTail = alphaTail;
point0.parameters.betaTail  = betaTail;
point0.parameters.etaAB     = etaAB;
point0.parameters.kAB       = kAB;
point0.parameters.fftM      = M;
point0.parameters.fftValue  = value;
point0.parameters.fftType   = type;
point0.parameters.xRange    = [-12, 12];
point0.parameters.nPoints   = 1000;

point0.pdf.x  = x;
point0.pdf.MA = fMA;
point0.pdf.GL = fGL;
point0.pdf.AB = fAB;

point0.pdfFigure  = fig;
point0.outputPath = '';

end