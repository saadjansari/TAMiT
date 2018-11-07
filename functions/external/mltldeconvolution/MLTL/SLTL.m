function SLTL(j, data, param, debug)

global MG;

% Step sizes for current level
tau = MG.pc.tau{j};

% Update the scaling-function coefficients at coarsest scale
if j == MG.J
	MG.es_hat{j+1} = tau(1) * MG.rs_hat{j+1};
	if param.lambdasf
		MG.es_hat{j+1} = MG.es_hat{j+1} - fftn(MG.Prox(param.lambdasf*tau(1), ifftn(MG.xs_hat{j+1} + MG.es_hat{j+1})));
	end
	MG.xs_hat{j+1} = MG.xs_hat{j+1} + MG.es_hat{j+1};
end

% Update the wavelet coefficients
MG.ew_hat{j} = cell(size(MG.rw_hat{j}));
for m = 1:numel(MG.ew_hat{j})
	MG.ew_hat{j}{m} = tau(m+1) * MG.rw_hat{j}{m};
	MG.ew_hat{j}{m} = MG.ew_hat{j}{m} - fftn(MG.Prox(param.lambda(j)*tau(m+1), ifftn(MG.xw_hat{j}{m} + MG.ew_hat{j}{m})));
	MG.xw_hat{j}{m} = MG.xw_hat{j}{m} + MG.ew_hat{j}{m};
end