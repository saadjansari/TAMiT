function UpdateLevel(j, data, param, debug)

global MG;

% Wavelet decomposition of the residual
[MG.rs_hat{j+1}, MG.rw_hat{j}] = WavAnalysisSep(MG.rs_hat{j}, MG.pc.hstr_hat(j, :), MG.pc.gstr_hat(j, :), param.decimation);

% Multilevel cycle
for mu = 1:param.MG.mu
	for eta = 1:param.MG.eta1
		UpdateResidual(j, param);
		SLTL(j, data, param, debug);
	end
	if j < MG.J
		if ~isempty(MG.ew_hat{j})
			if ~isempty(MG.es_hat{j+1})
				UpdateResidual(j, param);
			else
				for m = 1:numel(MG.ew_hat{j})
					MG.rs_hat{j+1} = MG.rs_hat{j+1} - conj(MG.pc.c_hat{j+1}{m+1}) .* MG.ew_hat{j}{m};
				end
			end
		end
		MG.es_hat{j+1} = MG.xs_hat{j+1};
		UpdateLevel(j+1, data, param, debug);
		MG.es_hat{j+1} = MG.xs_hat{j+1} - MG.es_hat{j+1};
	end
	for eta = 1:param.MG.eta2
		UpdateResidual(j, param);
		SLTL(j, data, param, debug);
	end
end

% Transfer all modifications to fine scale
MG.rs_hat{j+1} = [];
MG.rw_hat{j} = [];
MG.es_hat{j+1} = [];
MG.ew_hat{j} = [];
MG.xs_hat{j} = WavSynthesisSep(MG.xs_hat{j+1}, MG.xw_hat{j}, MG.pc.h_hat(j, :), MG.pc.g_hat(j, :), param.decimation);