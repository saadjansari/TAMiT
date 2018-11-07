% demo.m
%
% Demonstration script for the MLTL algorithm.
%
% (c) Cedric Vonesch, Biomedical Imaging Group, EPFL, 2009.04.16-2009.07.01

SetPath();

%% Generate data
% Original signal (Cameraman image)
debug.xorig_hat = fftn(double(imread('cameraman.pgm', 'pgm')));
showimage(1, real(ifftn(debug.xorig_hat)), 'Original');

% PSF (9-by-9 uniform blur)
data.h_hat = zeros(size(debug.xorig_hat));
data.h_hat(1:9, 1:9) = 1/81;
data.h_hat = circshift(data.h_hat, [-4 -4]);
data.h_hat = fftn(data.h_hat);

% Noise realization (white Gaussian with variance 0.308)
RandStream.setDefaultStream(RandStream('mt19937ar', 'seed', 0));
debug.sigma = sqrt(0.308);
b_hat = fftn(debug.sigma*randn(size(debug.xorig_hat)));

% Measurement
data.y_hat = data.h_hat .* debug.xorig_hat + b_hat;
showimage(2, real(ifftn(data.y_hat)), 'Blurred and noisy');

%% Algorithm parameters
% Wavelet family
param.wavelet = 'Haar';
% param.wavelet = 'Daub2';
% param.wavelet = 'Sym4';
% param.wavelet = '97';

% Number of decomposition levels (for each dimension)
param.J = [3 3];

% Decimated or undecimated (redundant) transform
param.decimation = true;

% Random shifts
% param.shift = 'none';
param.shift = 'random';
% param.shift = 'cycle';

% Regularization parameter for the coarsest-scale scaling-function subband
param.lambdasf = 0;

% Regularization parameters for the wavelet subbands (level-dependent)
param.lambda = [0.04 0.04 0.04];

% Thresholding function
param.threshold = 'Soft';
% param.threshold = 'Hard';

% Initial estimate
% param.iniest = 'zero';
% param.iniest = 'meas';
param.iniest = 'tik';

% Number of iterations
param.K = 50;

% Type of multilevel cycles
param.cycle = 'C'; % Coarse-to-fine
% param.cycle = 'V'; % V-cycles
% param.cycle = 'W'; % W-cycles
switch param.cycle
	case 'C'
		param.MG.eta1 = 0;
		param.MG.mu = 1;
		param.MG.eta2 = 1;
	case 'V'
		param.MG.eta1 = 1;
		param.MG.mu = 1;
		param.MG.eta2 = 1;
	case 'W'
		param.MG.eta1 = 1;
		param.MG.mu = 2;
		param.MG.eta2 = 1;
end

%% Run algorithm and compute SNR improvement
MLTLPrepare(data, param, debug);
result = MLTL(data, param, debug);
showimage(3, real(ifftn(result.x_hat)), 'Deconvolved');
fprintf('SERG = %f\n', ComputeSERGain(debug.xorig_hat, data.y_hat, result.x_hat));