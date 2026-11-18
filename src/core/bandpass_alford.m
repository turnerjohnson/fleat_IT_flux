function Xf = bandpass_alford(X, dt_hr, fc_cpd, c)
% BANDPASS_ALFORD  Zero-phase Butterworth bandpass, Alford (2003) Appendix A2.
%
%   Xf = bandpass_alford(X, dt_hr, fc_cpd, c)
%
%   X       [nrow x ntime]  filtered ALONG ROWS (time is dim 2). Must be gap-free.
%   dt_hr   sampling interval in hours (1 for your hourly-averaged record)
%   fc_cpd  centre frequency in cycles per day.  M2 = 1.9323.  f = 2*sind(lat).
%   c       bandwidth parameter. Alford uses 1.25: quarter-power points sit at
%           [fc/c, fc*c]. Narrow enough to isolate the band, wide enough that
%           the filter does not ring.
%
%   4th-order Butterworth designed on the half-band, run forwards and backwards
%   with filtfilt so there is no phase distortion. filtfilt doubles the effective
%   order, which is why the design order here is 4/2 = 2.
%
%   NOTE ON ORDER OF OPERATIONS: the modal fit is linear in eta at each timestep,
%   so bandpassing eta and then fitting gives exactly the same modal amplitudes as
%   fitting raw eta and then bandpassing the amplitude time series. Because your
%   eta has gaps but the amplitude series does not, fitting first is usually the
%   better route. This function is here for when you want banded eta for plotting.

if nargin < 4, c = 1.25; end

fs   = 24/dt_hr;                 % sampling rate, cycles per day
fnyq = fs/2;
Wn   = [fc_cpd/c, fc_cpd*c]/fnyq;

if Wn(2) >= 1
    error('bandpass_alford:nyquist', ...
          'Upper corner %.3f cpd is at or above Nyquist %.3f cpd.', fc_cpd*c, fnyq);
end

[b,a] = butter(2, Wn, 'bandpass');   % order 2 -> effective order 4 after filtfilt

if any(~isfinite(X(:)))
    error('bandpass_alford:gaps', ...
          'Input has NaNs. Fill them (fillmissing) before filtering.');
end

Xf = filtfilt(b, a, X.').';          % filtfilt works down columns, so transpose
end