function [c_hat, tau] = AnalyzeFilter(c_hat, h_hat, g_hat, decimation)
% [c_hat, tau] = AnalyzeFilter(f2_hat, h_hat, g_hat, decimation)
%
% Compute the correction filters and the step sizes for a filter f_hat,
% as required by the MLTL algorithm.
%
% Notes:
% - all filters are specified in the frequency (DFT) domain (_hat suffix);
% - the output arguments are cell arrays corresponding to the scales/subbands.
%
% f2_hat:       Squared modulus of the filter to be analyzed.
% h_hat, g_hat: Synthesis scaling and wavelet filters.
% decimation:   Decimated transform (true or false).
%
% c_hat:        Correction filters.
% tau:          Step sizes.
%
% (c) Cedric Vonesch, Biomedical Imaging Group, EPFL, 2006.12.14-2009.04.17

% Decomposition depth and dimensionality
J = size(h_hat, 1);

% Memory allocation and initializations
c_hat = [{{c_hat}}; cell(J, 1)];
rho = cell(J, 1);
tau = cell(J, 1);

% Compute the crosstalk constants and the equivalent filters within one pyramidal decomposition
for j = 1:J
	DecimateDim = ~cellfun(@isempty, h_hat(j, :));
	M = 2^(sum(DecimateDim));
	c_hat{j+1} = cell(M, 1);
	rho{j} = zeros(M);
	tau{j} = zeros(M, 1);
	m_hat = cell(M, M);

	% Upper triangular part of the crosstalk matrix
	for m1 = 1:M
		m1b = dec2bin(m1-1, log2(M));
		for m2 = m1:M
			m2b = dec2bin(m2-1, log2(M));

			tmp_hat = c_hat{j}{1};

			db = sum(DecimateDim)+1;
			for d = 1:numel(DecimateDim)
				if DecimateDim(d)
					db = db - 1;
					if m1b(db) == '0' && m2b(db) == '0'
						tmp_hat = WavAnalysis1D(tmp_hat, abs(h_hat{j, d}).^2, d, decimation);
					elseif m1b(db) == '0' && m2b(db) == '1'
						tmp_hat = WavAnalysis1D(tmp_hat, h_hat{j, d} .* conj(g_hat{j, d}), d, decimation);
					elseif m1b(db) == '1' && m2b(db) == '0'
						tmp_hat = WavAnalysis1D(tmp_hat, g_hat{j, d} .* conj(h_hat{j, d}), d, decimation);
					elseif m1b(db) == '1' && m2b(db) == '1'
						tmp_hat = WavAnalysis1D(tmp_hat, abs(g_hat{j, d}).^2, d, decimation);
					end
				end
			end
			rho{j}(m1, m2) = max(abs(tmp_hat(:)));
			if m1 == 1
				c_hat{j+1}{m2} = tmp_hat;
			end
			m_hat{m1, m2} = tmp_hat;
		end
	end

	% Lower triangular part of the crosstalk matrix (completed by symmetry)
	for m1 = 2:M
		for m2 = 1:m1-1
			rho{j}(m1, m2) = rho{j}(m2, m1);
		end
	end
end

% Compute the bounds
for j = 1:J-1
	for m = 2:numel(tau{j})
		tau{j}(m) = 1/sum(rho{j}(m, 2:end));
	end
end
for m = 1:numel(tau{J})
	tau{J}(m) = 1/sum(rho{J}(m, :));
end