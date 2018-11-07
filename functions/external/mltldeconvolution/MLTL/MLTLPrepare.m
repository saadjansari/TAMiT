function MLTLPrepare(data, param, debug)
% MLTLPrepare(data, param, debug)
%
% Perform initialization tasks before each call to MLTL.
%
% NOTE: all signals are specified in the frequency (DFT) domain,
%       which is indicated with a _hat suffix.
%
% data.y_hat:           Measured signal.
% data.h_hat:           Filter modeling the imaging system.
% param.decimation:     Whether the wavelet representation is decimated (true or false).
% param.threshold:      Thresholding function ('Soft' or 'Hard').
% param.shift:          Optionally shift the signal before every mu-cycle. See GenerateShifts() routine.
% param.iniest:         Initial estimate. See InitialEstimate() routine.
%
% (c) Cedric Vonesch, Biomedical Imaging Group, EPFL, 2008.04.03-2009.04.16

% Multigrid heap
global MG;
MG.J = max(param.J);
MG.xs_hat = cell(MG.J+1, 1);
MG.xw_hat = cell(MG.J, 1);
MG.es_hat = cell(MG.J+1, 1);
MG.ew_hat = cell(MG.J, 1);
MG.rs_hat = cell(MG.J+1, 1);
MG.rw_hat = cell(MG.J, 1);

% Precompute wavelet filters, correction filters and step sizes
[MG.pc.h_hat, MG.pc.g_hat, MG.pc.htld_hat, MG.pc.gtld_hat] = BiorthWavFiltersSep(GetDimensions(data.y_hat), param.J, param.decimation, param.wavelet);
MG.pc.hstr_hat = cellfuncell(@conj, MG.pc.h_hat);
MG.pc.gstr_hat = cellfuncell(@conj, MG.pc.g_hat);
MG.pc.htldstr_hat = cellfuncell(@conj, MG.pc.htld_hat);
MG.pc.gtldstr_hat = cellfuncell(@conj, MG.pc.gtld_hat);
[MG.pc.c_hat, MG.pc.tau] = AnalyzeFilter(abs(data.h_hat).^2, MG.pc.h_hat, MG.pc.g_hat, param.decimation);

% Reblurred measurement
MG.hstry_hat = conj(data.h_hat) .* data.y_hat;

% Thresholding function
MG.Prox = str2func(sprintf('Prox%s', param.threshold));

% Shift sequence
[MG.vargrid, MG.n0] = GenerateShifts(param.shift, GetDimensions(data.y_hat), param.K);

% Initial estimate and its wavelet decomposition
MG.xs_hat{1} = InitialEstimate(data, param, debug);
if ~MG.vargrid
	for j = 1:MG.J
		[MG.xs_hat{j+1}, MG.xw_hat{j}] = WavAnalysisSep(MG.xs_hat{j}, MG.pc.htldstr_hat(j, :), MG.pc.gtldstr_hat(j, :), param.decimation);
	end
end