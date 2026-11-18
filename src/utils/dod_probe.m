function info = dod_probe(datadir)
% DOD_PROBE  Say what every field in the M1-M5 DOD files actually is.
%
%   info = dod_probe()
%   info = dod_probe(datadir)
%
%   Prints, for each mooring, the shape / units / range of every variable in
%   M*_TSPbar.mat and M*_TZdepth.mat and in FLEAT_M_ADCP_8m_BindDp.mat, the
%   sample interval of each clock, and the three clock-alignment facts that
%   decide how the hourly grid has to be built.
%
%   info  struct of the numbers the print-out is built from, so a driver can
%         assert on them instead of trusting a remembered value.
%
%   WHAT THE FIELDS ARE  (all five moorings, verified not assumed)
%
%   M*_TZdepth.mat
%     MTime    [1 x 44569]  datenum, 10.000 min exactly, no gaps, monotonic.
%                           2016-05-25 00:00 -> 2017-03-30 12:00 (309.5 d).
%                           Byte-identical in TZdepth and TSPbar, and identical
%                           across all five moorings.
%     Zdepth   [1 x 95]     60:10:1000 m. A FIXED grid -- the gridded product,
%                           not instrument depths.
%     TZdepth  [44569 x 95] TEMPERATURE, degC, on that grid. Note the shape:
%                           TIME DOWN ROWS, depth across columns -- transposed
%                           relative to the fleat.F* convention of [nz x nt].
%
%   M*_TSPbar.mat   (columns are the T/S instruments, shallow -> deep)
%     *** THE INSTRUMENT COUNT IS NOT THE SAME AT EVERY MOORING ***
%         M1 5, M2 6, M3 6, M4 5, M5 5. Nothing may assume 5.
%     PSdbar   [44569 x ns] pressure at each instrument, dbar. Time-varying:
%                           this is the knockdown record. M1 nominal depths
%                           ~62 / 89 / 129 / 259 / 763 dbar.
%     TSdbar   [44569 x ns] temperature AT THE S INSTRUMENTS, degC. Not the
%                           same thing as TZdepth -- point sensors, not a
%                           grid. Use it with Sdbar, not instead of TZdepth.
%     Sdbar    [44569 x ns] salinity, psu.
%     Cdbar    [44569 x ns] conductivity, mS/cm. Redundant with S and T; kept
%                           only for provenance, and it is the field with the
%                           worst bad block (see below).
%
%   KNOWN BAD DATA  (small, but it will poison a mean or an S(T) fit)
%     M2 Cdbar col 3   876 CONTIGUOUS samples (~6.1 d, index 13436-14311) of
%                      C = -25, i.e. a fill value, not a measurement. Sdbar
%                      over the same span looks normal, so S was evidently
%                      repaired independently -- use Sdbar, not Cdbar.
%     Sdbar everywhere 126-306 isolated high spikes per mooring (0.06-0.11%),
%                      reaching 37.05 psu, which this region does not do.
%                      Mostly the upper sensors. These are spikes, not signal.
%     M3 col 2         4.09% NaN across TSdbar/Sdbar/Cdbar -- a sensor gap.
%     TZdepth          clean: nothing outside [3, 32] degC at any mooring, and
%                      inversions > 0.05 degC occur in under 0.05% of pairs.
%
%   FLEAT_M_ADCP_8m_BindDp.mat   (all five moorings in one file, suffix = m)
%     ADCPtm<m>      [1 x ~7540]  datenum, 60.000 min exactly. ALREADY HOURLY.
%     ZZ             [1 x 75]     8:8:600 m. Shared by all five. The "8m Bin"
%                                 in the filename; "Dp" = binned on depth.
%     U<m>, V<m>     [75 x ~7540] velocity in CM/S, not m/s. |u| p99 ~ 31-44.
%     dnADCPdepth<m> [1 x ~7540]  depth of the down-looking head, m. Varies
%                                 50 -> 190 m: knockdown again.
%     upADCPtilt<m>  [1 x ~7540]  *** NOT TILT. This is INSTRUMENT TEMPERATURE
%     dnADCPtilt<m>  [1 x ~7540]  *** in degC. See below.
%     ans            [1 x 1]      8.008. A stray variable left in the save by
%                                 whoever wrote the file -- the bin size. It
%                                 carries no record and should be ignored.
%
%   WHY THE TILT FIELDS ARE TEMPERATURE
%   Four independent reasons, all checked here by check_tilt below:
%     1. The values sit at 15-30 with a median of 28.7. A mooring that held a
%        28 deg list for ten months would not be a mooring.
%     2. They are ANTI-correlated with the instrument's own depth (r = -0.80 at
%        M1). Real tilt goes UP with blowdown, not down. Temperature at a
%        pulled-down sensor goes down, which is the sign observed.
%     3. Interpolate TZdepth to the ADCP's own depth and time and the two
%        agree at r = +0.978 to +0.991 with a bias under 0.33 degC and an rms
%        of 0.2-0.6 degC, at every mooring that has enough overlap.
%     4. up and dn agree with each other at r = 0.98 and differ by ~0.1 degC,
%        which is two sensors in the same water, not two tilt axes.
%   So: do not use these as a knockdown or data-quality flag, and do not feed
%   them to anything expecting degrees. Use dnADCPdepth for knockdown.
%
%   IS THE ADCP FILE EVEN THESE MOORINGS?  YES, AND THE NUMBERING MATCHES
%   Worth asking, because DOD also deployed Barnys (B1-B3) on the reef top.
%   The ADCP file is the five DEEP moorings, and ADCP suffix m is the same
%   physical mooring as M<m>_TZdepth / M<m>_TSPbar. Evidence, strongest first:
%
%   1. KNOCKDOWN. If an ADCP and a T/S chain hang off the same wire they get
%      pulled down together. Correlate dnADCPdepth<m> against mooring n's
%      shallowest PSdbar, all 25 combinations (check_identity below):
%              n=1     n=2     n=3     n=4     n=5
%        m=1  0.995   0.078   0.427   0.134   0.028
%        m=2  0.080   0.998  -0.071   0.030   0.499
%        m=3  0.439  -0.048   0.988   0.448   0.027
%        m=4  0.137   0.034   0.417   0.971   0.351
%        m=5  0.032   0.495   0.026   0.360   0.991
%      Diagonal 0.97-0.998, best off-diagonal 0.42-0.50. Not ambiguous.
%
%   2. THE OFF-DIAGONAL RECOVERS THE ARRAY GEOMETRY. The leftover correlation
%      is between moorings the paper puts near each other: 1-3 (0.43, both
%      west of the reef in >2600 m), 2-5 (0.50, both east), and 4 with 3 and 5
%      (0.45, 0.35 -- M4 sits north, between M3 and M5). The off-diagonal is
%      shared forcing across the array, and its pattern is the map. So the
%      numbering is the paper's numbering, not just internally consistent.
%
%   3. The instrument temperature in the "tilt" fields matches mooring m's own
%      T record better than any other mooring's at m = 1,2,3,4 (r 0.978-0.991
%      on the diagonal, 0.76-0.94 off it). M5 cannot be tested this way --
%      its T field is ~100% NaN at 60 m, right where its ADCP sits.
%
%   4. It cannot be the Barnys. There are five records here and only three
%      Barnys; every mooring has an up AND a down field pair, while a Barny
%      carries one upward-looking ADCP; and ZZ reaches 600 m, which is below
%      the 53-135 m the Barnys sit in.
%
%   5. ZZ = 8:8:600 matches the deep moorings' stated design, "to measure
%      current velocities in the upper 600 m", and the up/dn pairs match the
%      stated 300 kHz up + 75 kHz Long Ranger down in one buoy.
%
%   One loose end: the paper puts the ADCP buoys at 64-91 m, and the 1st
%   percentile of dnADCPdepth gives ~50 (M1), 66 (M2), 58 (M3), 85 (M4), 60
%   (M5) m. M1, M3 and M5 read 4-14 m shallower than the quoted range. That is
%   a depth-reference question, not an identity question, given point 1.
%
%   THE THREE FACTS THAT CONSTRAIN THE HOURLY GRID
%   a) 44569 = 6*7428 + 1. A straight 6-sample block mean (cfg.avg_k = 6 in
%      mooring_hourly_mean) drops one sample and lands on :25 past the hour,
%      not on the hour.
%   b) The ADCP clocks are hourly but OFF the hour, by a DIFFERENT amount at
%      each mooring: M1 +34.8, M2 +15.7, M3 +14.2, M4 +55.0, M5 +35.0 minutes.
%      So the T/S clock and the velocity clock never coincide, and no block
%      mean can make them. One of the two has to be put on a common grid.
%   c) The ADCP records start 1-4 days BEFORE and end 1-4 days AFTER the T/S
%      record. The usable window is the intersection, which is the T/S span.
%
%   THERE IS NO COMMON VERTICAL GRID EITHER
%      T  60:10:1000 m (95 levels, fixed)
%      S  5 instruments on time-varying pressure, ~60 to ~760 dbar
%      u,v 8:8:600 m (75 bins, fixed)
%   T, S and velocity all exist only over 60-600 m. build_dod_mooring takes
%   that overlap as its default grid.

%% ---- arguments ------------------------------------------------------
if nargin < 1 || isempty(datadir)
    datadir = '/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/';
end
if ~isfolder(datadir)
    error('dod_probe:nodir', 'Not a folder: %s', datadir);
end
adcpfile = fullfile(datadir, 'FLEAT_M_ADCP_8m_BindDp.mat');

fprintf('\n============ DOD M1-M5: what is in these files ============\n');
fprintf('datadir: %s\n', datadir);

%% ---- the T/S files -------------------------------------------------
for m = 1:5
    fz = fullfile(datadir, sprintf('M%d_TZdepth.mat', m));
    fb = fullfile(datadir, sprintf('M%d_TSPbar.mat',  m));
    if ~isfile(fz) || ~isfile(fb)
        fprintf('\n-- M%d: files missing, skipped\n', m); continue
    end
    Z = load(fz);  B = load(fb);

    fprintf('\n---- M%d ----\n', m);
    describe('TZdepth (T, degC)', Z.TZdepth);
    describe('Zdepth  (m)',       Z.Zdepth);
    fprintf('    Zdepth = %g:%g:%g m, %d levels, spacing unique = %s\n', ...
        Z.Zdepth(1), Z.Zdepth(2)-Z.Zdepth(1), Z.Zdepth(end), ...
        numel(Z.Zdepth), mat2str(unique(diff(Z.Zdepth))));
    for f = {'PSdbar','TSdbar','Sdbar','Cdbar'}
        describe(f{1}, B.(f{1}));
        fprintf('      per-column median: %s\n', ...
            num2str(median(B.(f{1}), 1, 'omitnan'), '%10.3f'));
    end

    info.M(m).clock = clockfacts(sprintf('M%d MTime', m), Z.MTime);
    fprintf('    MTime identical in TZdepth and TSPbar: %d\n', ...
        isequal(Z.MTime, B.MTime));
    if m > 1 && ~isempty(info.M(1).clock)
        fprintf('    MTime identical to M1: %d\n', isequal(Z.MTime, info.M(1).clock.t));
    end
    info.M(m).nt     = numel(Z.MTime);
    info.M(m).nsens  = size(B.PSdbar, 2);
    info.M(m).Zdepth = Z.Zdepth;
    fprintf('    >>> %d T/S instruments at this mooring <<<\n', info.M(m).nsens);
    info.M(m).bad = badblocks(B, info.M(m).clock.t);
end

%% ---- the ADCP file -------------------------------------------------
if ~isfile(adcpfile)
    fprintf('\nADCP file missing: %s\n', adcpfile);
    return
end
A = load(adcpfile);
fprintf('\n---- ADCP (all five moorings, one file) ----\n');
describe('ZZ (m)', A.ZZ);
fprintf('    ZZ = %g:%g:%g m, %d bins, spacing unique = %s\n', ...
    A.ZZ(1), A.ZZ(2)-A.ZZ(1), A.ZZ(end), numel(A.ZZ), mat2str(unique(diff(A.ZZ))));
if isfield(A,'ans')
    fprintf('    stray variable "ans" = %g  (bin size; ignore it)\n', A.ans);
end
for m = 1:5
    if ~isfield(A, sprintf('U%d', m)), continue; end
    fprintf('  -- mooring %d\n', m);
    describe(sprintf('U%d (cm/s)', m), A.(sprintf('U%d', m)));
    describe(sprintf('V%d (cm/s)', m), A.(sprintf('V%d', m)));
    describe(sprintf('dnADCPdepth%d (m)', m), A.(sprintf('dnADCPdepth%d', m)));
    describe(sprintf('upADCPtilt%d (degC!)', m), A.(sprintf('upADCPtilt%d', m)));
    describe(sprintf('dnADCPtilt%d (degC!)', m), A.(sprintf('dnADCPtilt%d', m)));
    info.ADCP(m).clock = clockfacts(sprintf('ADCPtm%d', m), A.(sprintf('ADCPtm%d', m)));
    % per-bin coverage: which bins are real
    U  = A.(sprintf('U%d', m));
    nf = mean(isnan(U), 2);
    gb = find(nf < 0.5);
    fprintf('      bins >50%% good: ZZ %g -> %g m (%d of %d). worst bin %.0f%% NaN\n', ...
        A.ZZ(gb(1)), A.ZZ(gb(end)), numel(gb), numel(A.ZZ), 100*max(nf));
    info.ADCP(m).nanfrac_by_bin = nf;
end

%% ---- clock alignment, the part that matters ------------------------
fprintf('\n---- clock alignment ----\n');
for m = 1:5
    if ~isfield(info,'M') || numel(info.M) < m || isempty(info.M(m).nt), continue; end
    tm = info.M(m).clock.t;  ta = info.ADCP(m).clock.t;
    offmin = mod((ta(1) - floor(ta(1)))*1440, 60);
    fprintf(['M%d: T/S dt=%.4f min, n=%d, n/6=%.4f (block mean drops %d, lands on ' ...
             ':%02.0f past)\n'], m, info.M(m).clock.dt_min, numel(tm), ...
             numel(tm)/6, numel(tm) - 6*floor(numel(tm)/6), ...
             mod((tm(1) + 2.5/1440 - floor(tm(1)))*1440, 60));
    fprintf('    ADCP dt=%.4f min, n=%d, sits %.2f min past the hour\n', ...
        info.ADCP(m).clock.dt_min, numel(ta), offmin);
    fprintf('    ADCP leads T/S start by %+.2f d, trails its end by %+.2f d\n', ...
        tm(1) - ta(1), ta(end) - tm(end));
    fprintf('    overlap %s -> %s (%.2f d)\n', ...
        datestr(max(tm(1), ta(1))), datestr(min(tm(end), ta(end))), ...
        min(tm(end), ta(end)) - max(tm(1), ta(1)));
    info.overlap(m,:) = [max(tm(1), ta(1)) min(tm(end), ta(end))];
end

%% ---- the tilt-is-temperature check --------------------------------
fprintf('\n---- is *ADCPtilt* tilt, or temperature? ----\n');
info.tilt = check_tilt(datadir, A);

fprintf('\n---- is ADCP suffix m the same mooring as M<m>_TSPbar? ----\n');
info.identity = check_identity(datadir, A);

fprintf('\nSuggested common vertical grid: 60:10:600 m (T, S and u,v all exist)\n');
fprintf('Next step: build_dod_mooring(m) -- hourly, on an exact-hour clock.\n\n');
end

%% ======================================================================
function bad = badblocks(B, t)
% Range-check each instrument record and say whether what fails is a spike or
% a block. The distinction decides the fix: a spike gets dropped, a block that
% long means the sensor or the processing failed and the whole span is suspect.
lim = struct('Sdbar', [32 35.6], 'Cdbar', [2.5 6.5], ...
             'TSdbar', [3 32],    'PSdbar', [10 1200]);
bad = struct();
for f = fieldnames(lim).'
    v = double(B.(f{1}));
    r = lim.(f{1});
    m = isfinite(v) & (v < r(1) | v > r(2));
    bad.(f{1}) = m;
    if ~any(m(:)), continue; end
    fprintf('    BAD %s: %d of %d (%.4f%%) outside [%g %g]\n', ...
        f{1}, sum(m(:)), numel(v), 100*mean(m(:)), r(1), r(2));
    for c = find(any(m, 1))
        i = find(m(:,c));
        fprintf('      col %d: n=%d  values [%.3f %.3f]  ', ...
            c, numel(i), min(v(i,c)), max(v(i,c)));
        if numel(i) == i(end) - i(1) + 1
            fprintf('CONTIGUOUS BLOCK idx %d-%d, %s -> %s (%.2f d)\n', ...
                i(1), i(end), datestr(t(i(1))), datestr(t(i(end))), ...
                t(i(end)) - t(i(1)));
        else
            fprintf('scattered spikes\n');
        end
    end
end
% a NaN gap big enough to matter
for f = {'TSdbar','Sdbar'}
    v = double(B.(f{1}));
    nf = mean(isnan(v), 1);
    for c = find(nf > 0.01)
        fprintf('    GAP %s col %d: %.2f%% NaN\n', f{1}, c, 100*nf(c));
    end
end
end

%% ======================================================================
function describe(name, v)
% One line: shape, NaN fraction, range. The units are in the caller's label
% because they are not recorded anywhere in the files.
v = double(v);
g = v(isfinite(v));
if isempty(g)
    fprintf('    %-26s %-14s ALL NaN/empty\n', name, mat2str(size(v)));
else
    fprintf('    %-26s %-14s NaN %5.1f%%  [%9.4g %9.4g] med %9.4g\n', ...
        name, mat2str(size(v)), 100*mean(~isfinite(v(:))), min(g), max(g), median(g));
end
end

%% ======================================================================
function c = clockfacts(name, t)
% Sample interval and whether the record is clean. Reports the MODE of dt as
% well as the median, because one odd step in a quarter-million is exactly the
% thing a median hides and a reshape-based block mean trips over.
t  = double(t(:)).';
dt = diff(t);
dt = dt(isfinite(dt));
c.t          = t;
c.dt_min     = median(dt)*1440;
c.dt_min_mode= mode(round(dt*86400, 3))/60;
c.dt_min_max = max(dt)*1440;
c.dt_min_min = min(dt)*1440;
c.monotonic  = all(dt > 0);
c.n_offstep  = sum(abs(dt - median(dt)) > 1e-9);
c.per_hour   = 1/(median(dt)*24);
fprintf(['    %-14s n=%6d  %s -> %s  %.3f d  dt med %.4f / mode %.4f min  ' ...
         'off-step %d  monotonic %d  per hour %.4f\n'], ...
    name, numel(t), datestr(t(1)), datestr(t(end)), t(end)-t(1), ...
    c.dt_min, c.dt_min_mode, c.n_offstep, c.monotonic, c.per_hour);
end

%% ======================================================================
function R = check_tilt(datadir, A)
% Compare dnADCPtilt against TZdepth sampled at the ADCP's own depth and time.
% If the field is temperature the two match; if it is tilt they cannot.
R = struct('m', {}, 'n', {}, 'r', {}, 'bias', {}, 'rms', {}, 'r_depth', {});
for m = 1:5
    fz = fullfile(datadir, sprintf('M%d_TZdepth.mat', m));
    if ~isfile(fz) || ~isfield(A, sprintf('ADCPtm%d', m)), continue; end
    Z   = load(fz);
    ta  = A.(sprintf('ADCPtm%d',      m));
    dep = A.(sprintf('dnADCPdepth%d', m));
    tlt = A.(sprintf('dnADCPtilt%d',  m));

    % nearest T/S sample in time, then interp down the fixed depth grid
    it   = max(1, min(numel(Z.MTime), interp1(Z.MTime, 1:numel(Z.MTime), ta, 'nearest', 'extrap')));
    Tcol = nan(size(ta));
    for k = 1:numel(ta)
        if isfinite(dep(k))
            Tcol(k) = interp1(Z.Zdepth, Z.TZdepth(it(k), :), dep(k), 'linear', NaN);
        end
    end
    ok = isfinite(Tcol) & isfinite(tlt);
    rd = corr_(dep(isfinite(dep) & isfinite(tlt)), tlt(isfinite(dep) & isfinite(tlt)));
    if sum(ok) > 100
        d = tlt(ok) - Tcol(ok);
        R(end+1) = struct('m', m, 'n', sum(ok), 'r', corr_(Tcol(ok), tlt(ok)), ...
            'bias', mean(d), 'rms', std(d), 'r_depth', rd);  %#ok<AGROW>
        fprintf(['M%d: n=%5d  corr(tilt, T@its own depth) = %+.4f  ' ...
                 'bias %+.3f degC  rms %.3f degC   corr(tilt, depth) = %+.3f\n'], ...
            m, R(end).n, R(end).r, R(end).bias, R(end).rms, rd);
    else
        fprintf(['M%d: only %d overlapping points (the ADCP sits above the 60 m ' ...
                 'top of Zdepth most of the time)  corr(tilt, depth) = %+.3f\n'], ...
            m, sum(ok), rd);
    end
end
if ~isempty(R) && median([R.r]) > 0.9
    fprintf(['VERDICT: *ADCPtilt* is INSTRUMENT TEMPERATURE in degC, median ' ...
             'r = %.3f against independent T. It is not tilt.\n'], median([R.r]));
end
end

function R = check_identity(datadir, A)
% Cross-correlate every ADCP's knockdown against every T/S mooring's
% knockdown. Same wire => same blowdown, so the diagonal has to win. The
% off-diagonal is not noise: it is shared forcing between nearby moorings, and
% its pattern should reproduce the array's layout.
P = cell(1,5);
for n = 1:5
    fb = fullfile(datadir, sprintf('M%d_TSPbar.mat', n));
    if isfile(fb)
        B = load(fb, 'MTime', 'PSdbar');
        P{n} = {double(B.MTime(:)).', double(B.PSdbar(:,1)).'};
    end
end
R = nan(5,5);
for m = 1:5
    if ~isfield(A, sprintf('ADCPtm%d', m)), continue; end
    ta  = double(A.(sprintf('ADCPtm%d', m)));
    dep = double(A.(sprintf('dnADCPdepth%d', m)));
    g   = isfinite(dep) & dep > 20;        % the record has a few <20 m spikes
    dg  = dep(g);  tg = ta(g);
    for n = 1:5
        if isempty(P{n}), continue; end
        pi_ = interp1(P{n}{1}, P{n}{2}, tg, 'linear', NaN);
        ok  = isfinite(pi_);
        if sum(ok) > 100, R(m,n) = corr_(dg(ok), pi_(ok)); end
    end
end
fprintf('        ');  fprintf('    n=%d ', 1:5);  fprintf('\n');
for m = 1:5
    fprintf('  m=%d ', m);  fprintf('%+8.3f', R(m,:));
    [~, b] = max(R(m,:));
    if b == m
        fprintf('   diagonal wins by %+.3f\n', R(m,m) - max(R(m,[1:m-1 m+1:5])));
    else
        fprintf('   *** best is n=%d, NOT the diagonal -- look ***\n', b);
    end
end
if all(diag(R)' == max(R, [], 2)')
    fprintf(['VERDICT: ADCP suffix m is the same mooring as M<m>_TSPbar at all ' ...
             'five. The off-diagonal is the array geometry.\n']);
end
end

%% ======================================================================
function r = corr_(a, b)
a = a(:); b = b(:);
r = NaN;
if numel(a) > 2
    c = cov(a, b);  r = c(1,2)/sqrt(c(1,1)*c(2,2));
end
end
