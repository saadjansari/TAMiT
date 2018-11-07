function SERG = ComputeSERGain(xorig, y, x)
% SERG = ComputeSERGain(xorig, y, x)
%
% Compute a signal-to-error-ratio gain.
%
% xorig: Original signal.
% y:     Measurement.
% x:     Estimate of the original signal.
%
% SERG:  Signal-to-error ratio gain in dB.
%
% (c) Cedric Vonesch, Biomedical Imaging Group, EPFL, 2007.02.11-2009.04.17

SERG = 20 * log10(norm(xorig(:)-y(:)) / norm(xorig(:)-x(:)));