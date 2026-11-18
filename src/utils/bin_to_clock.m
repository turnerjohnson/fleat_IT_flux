function [Y, n] = bin_to_clock(t, X, tc)
% BIN_TO_CLOCK  Centred bin mean of X onto an evenly spaced target clock.
%
%   [Y, n] = bin_to_clock(t, X, tc)
%
%   t    [1 x nt]        source timestamps, datenum
%   X    [nrow x nt]     data on that clock ([nt x nrow] is tolerated and
%                        transposed, same as build_dod_mooring's local hbin)
%   tc   [1 x nc]        target clock, evenly spaced. NOT assumed to be exact
%                        hours and NOT assumed to start on the hour -- the
%                        FLEAT luecke ADCP clock is hourly but may sit off the
%                        hour, and the point of this function is to bin onto
%                        THAT clock rather than onto a clock of our own.
%
%   Y    [nrow x nc]     bin mean, NaNs omitted. A bin with no finite sample
%                        comes back NaN rather than 0.
%   n    [nrow x nc]     how many finite samples went into each value, so a
%                        partial end bin or a gap is visible as a count and
%                        does not have to be inferred from the value.
%
%   WHY NEAREST-TIMESTAMP AND NOT A LEADING BLOCK MEAN
%   Each sample goes to the nearest tc, so the bin around tc(k) is
%   [tc(k)-dt/2, tc(k)+dt/2) -- symmetric, hence a zero-phase low-pass. A
%   leading block mean (samples 1:6, 7:12, ...) is what reshape gives you and
%   it lands dt*5/12 late, which for 10-minute data in hourly bins is 25
%   minutes, i.e. 12 degrees of M2 phase thrown away for free. For internal-
%   tide flux a phase error rotates u'p' into the wrong quadrant, so it is the
%   expensive kind of error. See WHICH WAY TO MOVE THE ADCP in
%   build_dod_mooring for the same argument made the other direction.
%
%   WHY THIS IS NOT build_dod_mooring's hbin
%   That one hard-codes exact UTC hours (idx = round((t-t_h(1))*24)+1) because
%   at the DOD moorings it builds the clock itself and can make it exact. Here
%   the clock is handed to us by somebody else's file, so the step comes from
%   the target vector. build_dod_mooring is deliberately left alone: its clock
%   is exact by construction, so routing it through here would buy nothing and
%   would touch a working pipeline file in a data-prep pass.
%
%   author TEJ, with help from Claude -- 17 November 2026

t  = double(t(:)).';
tc = double(tc(:)).';
nt = numel(t);
nc = numel(tc);
if nc < 2
    error('bin_to_clock:shortclock', ...
        'The target clock needs at least 2 samples to have a step; got %d.', nc);
end
if size(X, 2) ~= nt
    if size(X, 1) == nt
        X = X.';                                % tolerate [nt x nrow]
    else
        error('bin_to_clock:shape', ...
            ['X is %s and the source clock has %d samples, so neither ' ...
             'dimension of X is time. Nothing here can guess which it is.'], ...
             mat2str(size(X)), nt);
    end
end

dt  = median(diff(tc));
off = abs(diff(tc) - dt);
if any(off > 1e-6)                              % 1e-6 d = 0.086 s
    error('bin_to_clock:uneven', ...
        ['The target clock is not evenly spaced (%d steps differ from the ' ...
         'median %.6f d by more than 0.09 s, worst %.3f s). A nearest-index ' ...
         'bin assignment is wrong on an uneven clock -- interpolate in time ' ...
         'instead, or feed an even clock.'], nnz(off > 1e-6), dt, max(off)*86400);
end

idx = round((t - tc(1))/dt) + 1;                % nearest target sample
ok  = idx >= 1 & idx <= nc;
M   = sparse(find(ok), idx(ok), 1, nt, nc);

G   = isfinite(X);
X0  = X;  X0(~G) = 0;
num = X0 * M;
n   = double(G) * M;
Y   = num ./ n;                                 % 0/0 -> NaN, which is right
Y(n == 0) = NaN;
end
