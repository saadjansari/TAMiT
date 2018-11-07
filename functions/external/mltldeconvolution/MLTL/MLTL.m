function result = MLTL(data, param, debug)
% result = MLTL(data, param, debug)
%
% Multilevel thresholded Landweber algorithm for multidimensional
% signal restoration with a wavelet-domain sparsity constraint.
%
% Note: run MLTLPrepare() before calling this routine!
%
% Reference:
% Vonesch, Unser. A fast multilevel algorithm for wavelet-regularized image
% restoration. IEEE Transactions on Image Processing, March 2009, 18(3), 509-523.
%
% INPUT ARGUMENTS:
% param.decimation:     Whether the wavelet representation is decimated (true or false).
% param.lambdasf:       Regularization parameter for the coarsest-scale scaling-function
%                       coefficients (set to zero if no thresholding is desired).
% param.lambda:         Level-dependent regularization parameter for the wavelet coefficients.
% param.K:              Number of mu-cycles.
% param.MG.eta1:        Number of pre-relaxations.
% param.MG.mu:          Number of recursive calls in UpdateLevel (V-cycle for mu=1; W-cycle for mu=2).
% param.MG.eta2:        Number of post-relaxations.
% debug.active:         Debug mode (true or false; see MLTLDebug() routine).
%
% OUTPUT ARGUMENTS:
% result.x_hat:         Estimate of the minimizer.
%
% (c) Cedric Vonesch, Biomedical Imaging Group, EPFL, 2006.11.15-2009.04.17

global MG;

% Main loop
for k = 1:param.K
	if MG.vargrid
		MG.xs_hat{1} = Modulate(MG.n0(:, k+1)-MG.n0(:, k), MG.xs_hat{1});
		for j = 1:MG.J
			[MG.xs_hat{j+1}, MG.xw_hat{j}] = WavAnalysisSep(MG.xs_hat{j}, MG.pc.htldstr_hat(j, :), MG.pc.gtldstr_hat(j, :), param.decimation);
		end
		MG.rs_hat{1} = Modulate(MG.n0(:, k+1), MG.hstry_hat) - MG.pc.c_hat{1}{1} .* MG.xs_hat{1};
	else
		MG.rs_hat{1} = MG.hstry_hat - MG.pc.c_hat{1}{1} .* MG.xs_hat{1};
	end
	UpdateLevel(1, data, param, debug);
end

% Clean up
if MG.vargrid
	result.x_hat = Modulate(-MG.n0(:, k+1), MG.xs_hat{1});
else
	result.x_hat = MG.xs_hat{1};	
end
clear('global', 'MG');