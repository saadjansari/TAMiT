function x_hat = InitialEstimate(data, param, debug)
% x_hat = InitialEstimate(data, param, debug)
%
% Return the initial estimate for an iterative algorithm,
% depending on the value of param.iniest:
% - 'zero': the null signal.
% - 'meas': the measured signal.
% - 'tik': Tikhonov estimate.
%
% Note: depending on the option, additional arguments may be required.
%
% (c) Cedric Vonesch, Biomedical Imaging Group, EPFL, 2008.04.03-2009.04.16

switch lower(param.iniest)
	case 'zero'
		x_hat = zeros(size(data.y_hat));
	case 'meas'
		x_hat = data.y_hat;
	case 'tik'
		gamma = 10^(-3) * debug.sigma^2;
		x_hat = conj(data.h_hat) .* data.y_hat ./ (abs(data.h_hat).^2 + gamma);
end