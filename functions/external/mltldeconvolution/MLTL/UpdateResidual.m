function UpdateResidual(j, param)

global MG;

if ~isempty(MG.ew_hat{j}) % We need to go back to the next finer level
	MG.rs_hat{j+1} = [];
	MG.rw_hat{j} = [];
	MG.rs_hat{j} = MG.rs_hat{j} - MG.pc.c_hat{j}{1} .* WavSynthesisSep(MG.es_hat{j+1}, MG.ew_hat{j}, MG.pc.h_hat(j, :), MG.pc.g_hat(j, :), param.decimation);
	MG.es_hat{j+1} = [];
	MG.ew_hat{j} = [];
	[MG.rs_hat{j+1}, MG.rw_hat{j}] = WavAnalysisSep(MG.rs_hat{j}, MG.pc.hstr_hat(j, :), MG.pc.gstr_hat(j, :), param.decimation);
elseif ~isempty(MG.es_hat{j+1}) % We can use intra-level correction
	for m = 1:numel(MG.rw_hat{j})
		MG.rw_hat{j}{m} = MG.rw_hat{j}{m} - MG.pc.c_hat{j+1}{m+1} .* MG.es_hat{j+1};
	end
end