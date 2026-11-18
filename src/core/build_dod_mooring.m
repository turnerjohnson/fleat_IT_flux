function [F, info] = build_dod_mooring(m, opt)
% BUILD_DOD_MOORING  One DOD mooring (M1-M5) as a fleat-style hourly struct.
%
%   F        = build_dod_mooring(m)
%   [F, info]= build_dod_mooring(m, opt)
%
%   m    mooring number 1-5.
%   opt  .datadir   where the M*_*.mat and FLEAT_M_ADCP_8m_BindDp.mat live
%        .zgrid     output depth grid, m. Default 10:10:1000 -- the UNION of
%                   the three ranges, so each variable is NaN outside its own.
%        .s_extend_deep  true (default). Hold S constant below the deepest
%                   instrument so N2 survives to 1000 m -- see below.
%        .s_mask_gap  true (default). With s_method 'ts', leave S NaN for
%                   water inside the ~6.9-8.0 degC void between instruments.
%        .uv_method 'interp' (default) | 'binmean'  how velocity gets onto the
%                   hourly clock. Read WHICH WAY TO MOVE THE ADCP before
%                   changing this; it is a phase-versus-amplitude choice.
%        .s_method  'pressure' (default) | 'ts'     how 5 salinity sensors
%                   become a gridded S field. Both are approximations.
%        .uv_units  'm/s' (default) | 'cm/s'        the files are cm/s.
%        .qc        true (default). Range-reject the known bad instrument
%                   values BEFORE averaging -- see KNOWN BAD DATA below.
%        .lat .lon  needed for N2 only. Default 7.9 N / 134.6 E, the FLEAT
%                   array off Velasco Reef -- PLACEHOLDERS, see N2 below.
%
%   F    struct on one exact-hour clock, shaped like fleat.F*:
%          time   [1 x nh]      datenum, exact UTC hours
%          z      [1 x nz]      the output grid
%          id     [1 x nz]      level labels (= z; these are grid levels, not
%                               instruments -- the product is already gridded)
%          depth  [nz x nh]     z replicated. Constant in time by construction
%          T      [nz x nh]     degC
%          S      [nz x nh]     psu
%          u, v   [nz x nh]     m/s (converted from the files' cm/s)
%          spd    [nz x nh]     hypot(u,v)
%          dTdz   [nz x nh]     degC/m, positive up
%          N2     [nz x nh]     s^-2, or all-NaN if gsw is not on the path
%          inst   struct        the T/S instruments, hourly, on their own
%                               time-varying pressure: .P .T .S .C [ns x nh].
%                               ns is 5 at M1/M4/M5 and 6 at M2/M3 -- it is
%                               read from the file, never assumed
%          adcp   struct        .depth [1 x nh] down-head depth (knockdown)
%                               .Tup .Tdn [1 x nh] instrument TEMPERATURE
%          readme cellstr       every choice made below, in the struct
%
%   info .nh .t_step_min .avg_k .adcp_offset_min .dropped and the coverage
%        fractions, so a driver can assert instead of trusting this header.
%
%   WHY THIS FUNCTION EXISTS RATHER THAN A CALL TO MOORING_HOURLY_MEAN
%   mooring_hourly_mean block-averages every field that is on F.time, which is
%   exactly right once there IS an F.time. Here there are three clocks and
%   three vertical grids and no F, so the job is to build one. Specifically:
%     - T/S are 10-minute, 44569 samples, starting exactly 2016-05-25 00:00.
%       44569 = 6*7428 + 1, so a 6-sample block mean drops a sample AND lands
%       on :25 past the hour.
%     - the ADCP is already hourly but sits 14 to 55 minutes past the hour, a
%       different offset at every mooring, so no block mean can align the two.
%     - T is on 60:10:1000 m, velocity on 8:8:600 m, S on 5-6 instruments
%       whose pressure moves with the knockdown.
%   Once F comes back from here it is a single-clock struct, and
%   mooring_hourly_mean / trim_mooring_record / psectime_mooring apply to it
%   normally. Run dod_probe for the evidence behind every number above.
%
%   THE CLOCK
%   Exact UTC hours, t_h = tm(1) : 1/24 : tm(end). tm(1) is already midnight
%   so no rounding is needed at the start. T/S is averaged in CENTRED hourly
%   bins -- the six samples nearest each hour, i.e. [:-30, :+30) -- which is a
%   symmetric low-pass with no phase shift, unlike a leading block mean.
%
%   WHICH WAY TO MOVE THE ADCP  (opt.uv_method)
%   Hourly data that needs to move by up to half an hour. Neither option is
%   free, and for internal-tide work the cost is not symmetric:
%     'interp'  linear in time. No phase error. Damps M2 by at most
%               1 - cos(pi*1h/12.42h) = 3%, worst at the half-hour offset.
%     'binmean' nearest hour. No amplitude loss, but a fixed phase error of up
%               to 360*0.5/12.42 = 14.5 deg of M2 -- which at M4 (55 min off
%               the hour) is near its worst.
%   A 14 deg phase error rotates u'p' into the wrong quadrant and biases the
%   tidal flux; a 3% amplitude loss scales it. So 'interp' is the default.
%   NaN handling differs too: interp1 spreads each NaN to its two neighbouring
%   hours, which matters at M4, where the upper 12 bins are 56% missing.
%
%   KNOWN BAD DATA, AND WHY IT IS REJECTED BEFORE THE MEAN  (opt.qc)
%   Three defects, all small, all able to do damage out of proportion to their
%   size because a 37 psu spike survives an hourly mean as a sixth of itself
%   and a -25 conductivity does not survive it at all:
%     - M2 Cdbar col 3 has 876 CONTIGUOUS samples (~6.1 d) of C = -25. That is
%       a fill value. Sdbar over the same span is normal, so S was repaired
%       independently and C was not. This is why C is provenance only here and
%       S is taken from Sdbar.
%     - every mooring has 126-306 isolated Sdbar high spikes (0.06-0.11% of
%       samples), up to 37.05 psu, mostly in the upper sensors. Left in, they
%       bend the S(T) curve of opt.s_method 'ts' and put spurious inversions
%       into N2.
%     - M3 col 2 is 4.09% NaN across T/S/C: a sensor gap, not a spike. It is
%       left as NaN and the pressure interpolation simply spans it.
%   So opt.qc NaNs anything outside S [32 35.6], C [2.5 6.5], T [3 32],
%   P [10 1200] BEFORE binning, and reports what it dropped. TZdepth itself is
%   clean -- nothing outside [3 32] degC anywhere, inversions > 0.05 degC in
%   under 0.05% of pairs -- so the same test there is a no-op kept for safety.
%
%   HOW SPARSE SALINITY SENSORS BECOME A FIELD  (opt.s_method)
%     'pressure' interpolate S across the instruments in pressure, hour by
%                hour, using their own PSdbar so the knockdown is respected.
%                Honest about what was measured; crude in the vertical, and it
%                cannot resolve anything the sensors do not straddle.
%     'ts'       fit S as a function of T from the pooled instrument records,
%                then map that curve through the 95-level T field. Gives S at
%                the resolution of T, and in this part of the western Pacific
%                the T-S relation is tight enough that this is the standard
%                move. It is still a relation, not a measurement: it cannot
%                see salinity variability at fixed T, which is precisely the
%                barrier-layer / fresh-lens signal the region is known for.
%   Above the shallowest sensor the 'pressure' method clamps to that sensor's
%   value (the default grid starts at 60 m, 2 dbar above it). Below the
%   deepest it returns NaN rather than extrapolating.
%
%   THE GRID IS THE UNION OF THE THREE RANGES, NOT THEIR INTERSECTION
%   Default 10:10:1000. Each variable carries its own true vertical extent and
%   is NaN outside it, because the deepest useful depth is different for every
%   one of them:
%       u, v   8 to 600 m.   A HARD limit: ZZ = 8:8:600 is the only depth
%                            vector in the ADCP file, and the deep bins are
%                            0.0% NaN at all five moorings, so 600 m is where
%                            the 75 kHz down-looker stops, not a quality cut.
%                            There is no velocity below 600 m to recover.
%       T      to 1000 m.    0.00% NaN at 1000 m at all five -- 400 m deeper
%                            than velocity.
%       S      to ~676 (M5) - 776 (M2) dbar, the deepest instrument.
%   An intersection grid (the 120:10:600 this used to default to) is right
%   only if every output column must be complete in T, S, u and v at once. It
%   is wrong for a modal fit that builds rho' and u' from separate, non-
%   colocated fits: it would throw away 400 m of good T for nothing.
%
%   NEAR THE SURFACE EVERYTHING IS NaN, AND THAT IS THE HONEST ANSWER
%   TZdepth NaN fraction by level -- the nominal 60 m top of Zdepth is not the
%   top of its data, because the knockdown carries the shallowest sensor below
%   the grid top for much of the record and M4's top sensor starts at 94 dbar:
%            60m    70m    80m    90m   100m   110m   120m
%      M1   39.5%  23.2%  14.1%   8.3%   4.4%   2.0%   0.9%
%      M2  100.0%  68.3%  34.8%  18.2%  10.2%   5.1%   2.0%
%      M3   58.8%  20.7%   7.9%   4.0%   3.3%   2.5%   1.7%
%      M4  100.0% 100.0% 100.0%  34.4%   6.9%   0.9%   0.0%
%      M5   99.8%  44.8%  21.9%  10.5%   3.8%   0.6%   0.1%
%   So T and S are NaN over roughly the top 90-120 m and nothing here invents
%   values for them. Velocity is NOT -- it runs to 8 m, which is the point: on
%   this grid the levels above ~100 m carry u and v with no T or S beside them.
%   One exception worth knowing about, since it bites exactly there: M4's
%   upper 12 bins (8-96 m) are 55.7% NaN -- its up-looker partially failed.
%
%   S BELOW THE DEEPEST SENSOR IS HELD CONSTANT  (opt.s_extend_deep)
%   Default true, and it is not the same kind of choice as extrapolating S
%   upward. If S stops at its deepest instrument then N2 stops there too, and
%   220-320 m of perfectly good T at the bottom of the profile is lost -- the
%   bottom gappiness that hurts a mode fit most. The deep ocean makes this
%   cheap: S std at the deepest sensor is 0.0039 (M1), 0.0077 (M2), 0.0054
%   (M3), 0.0112 (M4), 0.0101 (M5) psu. Holding it constant costs ~0.01 psu,
%   about 0.008 kg/m3 of density, which is negligible against T's share of
%   deep N2. Set false to get NaN instead and see the truncation.
%
%   WHAT STILL IS NOT SOLVED BY ANY GRID CHOICE
%   Bottom depths at these moorings run 1350-3610 m, so 1000 m is 28-74% of
%   the water column. A modal decomposition over the full depth needs N2 below
%   the instrumented range from somewhere else -- climatology, argo_n2_monthly,
%   set_N2_source. Widening this grid cannot manufacture it.
%
%   N2 IS THE WEAKEST FIELD HERE AND THE POSITION IS A PLACEHOLDER
%   N2 needs gsw and a latitude. The DOD mooring positions are not in these
%   files and not in fleat_moorings, which carries F1-F6 and not M1-M5, so
%   opt.lat/opt.lon default to the FLEAT array position off Velasco Reef and
%   are almost certainly a few tenths of a degree wrong. That costs very
%   little in N2 (gravity and the SA conversion), but it is a guess, it is
%   recorded in F.readme, and it should be replaced from the deployment table
%   before N2 is used quantitatively. N2 also inherits whichever S
%   approximation was chosen, so at 'pressure' its vertical structure is only
%   as good as 5 or 6 points.

%% ---- arguments ------------------------------------------------------
if nargin < 1 || ~isscalar(m) || ~ismember(m, 1:5)
    error('build_dod_mooring:m', 'Needs a mooring number 1-5; got %s.', mat2str(m));
end
if nargin < 2 || ~isstruct(opt), opt = struct(); end
opt = setdef(opt, 'datadir',   '/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/');
opt = setdef(opt, 'zgrid',     10:10:1000);
opt = setdef(opt, 'uv_method', 'interp');
opt = setdef(opt, 's_method',  'pressure');
opt = setdef(opt, 'uv_units',  'm/s');
opt = setdef(opt, 'qc',        true);
opt = setdef(opt, 's_extend_deep', true);
opt = setdef(opt, 's_mask_gap',    true);
opt = setdef(opt, 'lat',       7.9);
opt = setdef(opt, 'lon',       134.6);

fz = fullfile(opt.datadir, sprintf('M%d_TZdepth.mat', m));
fb = fullfile(opt.datadir, sprintf('M%d_TSPbar.mat',  m));
fa = fullfile(opt.datadir, 'FLEAT_M_ADCP_8m_BindDp.mat');
for f = {fz, fb, fa}
    if ~isfile(f{1}), error('build_dod_mooring:nofile', 'Missing: %s', f{1}); end
end
Z = load(fz);  B = load(fb);
A = load(fa, sprintf('ADCPtm%d', m), sprintf('U%d', m), sprintf('V%d', m), ...
              sprintf('dnADCPdepth%d', m), sprintf('upADCPtilt%d', m), ...
              sprintf('dnADCPtilt%d', m), 'ZZ');

%% ---- check the clocks before using them ----------------------------
tm = double(Z.MTime(:)).';
if ~isequal(Z.MTime, B.MTime)
    error('build_dod_mooring:clockmismatch', ...
        ['M%d TZdepth and TSPbar are on different clocks (max diff %.3g s). ' ...
         'They are byte-identical in the five files checked, so this is new ' ...
         'and the two have to be interpolated separately.'], ...
         m, max(abs(Z.MTime - B.MTime))*86400);
end
dtm = diff(tm);
if any(abs(dtm - median(dtm)) > 1e-9)
    warning('build_dod_mooring:unevenTS', ...
        ['M%d MTime is not evenly spaced (%d off-step samples). The binning ' ...
         'below is by timestamp so it still works, but check for gaps.'], ...
         m, sum(abs(dtm - median(dtm)) > 1e-9));
end
info.t_step_min = median(dtm)*1440;
info.avg_k      = round(1/(median(dtm)*24));     % 6 samples per hour

ta = double(A.(sprintf('ADCPtm%d', m))(:)).';
info.adcp_step_min  = median(diff(ta))*1440;
info.adcp_offset_min= mod((ta(1) - floor(ta(1)))*1440, 60);

%% ---- the output clock: exact UTC hours over the T/S record ---------
t0 = ceil(tm(1)*24)/24;              % tm(1) is already midnight; harmless
t1 = floor(tm(end)*24)/24;
t_h = t0 : 1/24 : t1;
nh  = numel(t_h);
info.nh       = nh;
info.dropped  = numel(tm) - info.avg_k*floor(numel(tm)/info.avg_k);
z   = double(opt.zgrid(:)).';
nz  = numel(z);

%% ---- reject the known bad values BEFORE any averaging --------------
% Order matters. A spike removed after the mean has already been smeared over
% the other five samples in its hour; removed before, it costs one sample.
info.nsens = size(B.PSdbar, 2);                 % 5 at M1/M4/M5, 6 at M2/M3
[B, Z.TZdepth, info.qc] = qc_reject(B, double(Z.TZdepth), opt, m);

%% ---- T: 10-min -> centred hourly bins, then onto z -----------------
% TZdepth is [nt x nz]: time down ROWS. Transpose first -- this is the single
% easiest thing to get wrong in these files, and it does not error, it just
% silently gives a 95-hour record with 44569 depth levels.
T_full = hbin(tm, double(Z.TZdepth).', t_h);             % [95 x nh]
zz_T   = double(Z.Zdepth(:)).';
F.T    = interp1(zz_T, T_full, z, 'linear', NaN);

%% ---- the T/S instruments, hourly on their own pressure ------------
P_i = hbin(tm, double(B.PSdbar).', t_h);                 % [ns x nh]
T_i = hbin(tm, double(B.TSdbar).', t_h);
S_i = hbin(tm, double(B.Sdbar ).', t_h);
C_i = hbin(tm, double(B.Cdbar ).', t_h);

switch lower(opt.s_method)
    case 'pressure',  F.S = s_from_pressure(P_i, S_i, z, opt.s_extend_deep);
    case 'ts',        F.S = s_from_ts(T_i, S_i, F.T, opt.s_mask_gap);
    otherwise, error('build_dod_mooring:smethod', ...
            'opt.s_method must be ''pressure'' or ''ts''; got ''%s''.', opt.s_method);
end

%% ---- velocity: hourly-but-offset -> exact hours, then onto z -------
U = double(A.(sprintf('U%d', m)));                        % [75 x nta], cm/s
V = double(A.(sprintf('V%d', m)));
switch lower(opt.uv_method)
    case 'interp'
        U_h = interp1(ta, U.', t_h, 'linear', NaN).';
        V_h = interp1(ta, V.', t_h, 'linear', NaN).';
    case 'binmean'
        U_h = hbin(ta, U, t_h);
        V_h = hbin(ta, V, t_h);
    otherwise, error('build_dod_mooring:uvmethod', ...
            'opt.uv_method must be ''interp'' or ''binmean''; got ''%s''.', opt.uv_method);
end
zz_U = double(A.ZZ(:)).';
F.u  = interp1(zz_U, U_h, z, 'linear', NaN);
F.v  = interp1(zz_U, V_h, z, 'linear', NaN);
if strcmpi(opt.uv_units, 'm/s')
    F.u = F.u/100;  F.v = F.v/100;                        % files are cm/s
end
F.spd = hypot(F.u, F.v);

%% ---- derived ------------------------------------------------------
F.dTdz = nan(nz, nh);
if nz > 1
    for k = 1:nh, F.dTdz(:,k) = gradient(F.T(:,k), z); end
    F.dTdz = -F.dTdz;        % z is depth-positive-down; report degC/m up
end
[F.N2, n2note] = n2_of(F.T, F.S, z, opt);

%% ---- assemble in the fleat field order ----------------------------
F.time  = t_h;
F.z     = z;
F.id    = z;
F.depth = repmat(z(:), 1, nh);
F.inst  = struct('P', P_i, 'T', T_i, 'S', S_i, 'C', C_i);
F.adcp  = struct('depth', hbin(ta, double(A.(sprintf('dnADCPdepth%d', m))), t_h), ...
                 'Tup',   hbin(ta, double(A.(sprintf('upADCPtilt%d',  m))), t_h), ...
                 'Tdn',   hbin(ta, double(A.(sprintf('dnADCPtilt%d',  m))), t_h));
F = orderfields(F, {'time','z','id','depth','T','S','u','v','spd', ...
                    'dTdz','N2','inst','adcp'});

F.readme = { ...
 sprintf('M%d, built by build_dod_mooring on %s', m, datestr(now, 31))
 sprintf('time: %d exact UTC hours, %s to %s', nh, datestr(t_h(1)), datestr(t_h(end)))
 sprintf('T/S native step %.4f min (%d per hour); velocity native %.4f min, %.1f min past the hour', ...
         info.t_step_min, info.avg_k, info.adcp_step_min, info.adcp_offset_min)
 'T, S and the instrument records: centred hourly bin mean, omitnan, no phase shift'
 sprintf('u,v onto the hourly clock by ''%s''; see WHICH WAY TO MOVE THE ADCP in the header', opt.uv_method)
 sprintf('u,v in %s (the files hold cm/s)', opt.uv_units)
 sprintf('z = %g:%g:%g m -- the UNION of the three ranges, not their intersection', z(1), z(2)-z(1), z(end))
 'u,v exist 8-600 m (ADCP hard limit); T to 1000 m; S to the deepest instrument (~676-776 dbar)'
 'T and S are NaN over roughly the top 90-120 m: no sensor reaches there. u,v are not -- that is deliberate'
 sprintf('S below the deepest instrument: %s', deep_s_note(opt.s_extend_deep))
 'T from TZdepth (60:10:1000 m), interpolated to z. TZdepth in the file is [nt x nz] and is transposed here'
 sprintf('S by ''%s'' from %d instruments -- an approximation either way, see header', ...
         opt.s_method, info.nsens)
 sprintf('instrument count at this mooring: %d (5 at M1/M4/M5, 6 at M2/M3 -- read, not assumed)', ...
         info.nsens)
 info.qc.note
 sprintf('qc dropped: S %d, C %d, T %d, P %d samples', info.qc.dropped.Sdbar, ...
         info.qc.dropped.Cdbar, info.qc.dropped.TSdbar, info.qc.dropped.PSdbar)
 'M2 Cdbar col 3 holds ~6.1 d of C = -25 (fill value); Sdbar there is normal. C is provenance only'
 'depth is z replicated: the product is already gridded, so it does not move with the knockdown'
 'id is the level label (= z), NOT an instrument id -- these are grid levels'
 'adcp.Tup / adcp.Tdn come from the file''s upADCPtilt / dnADCPtilt fields, which hold'
 '  INSTRUMENT TEMPERATURE in degC, not tilt (r = +0.98 to +0.99 against TZdepth at'
 '  the instrument''s own depth, bias < 0.33 degC). Do not use them as a tilt or'
 '  knockdown flag -- adcp.depth is the knockdown record.'
 'the file''s stray variable "ans" (= 8.008, the bin size) is not loaded'
 n2note
 sprintf('N2 position used: %.4f N, %.4f E -- PLACEHOLDER, replace from the deployment table', ...
         opt.lat, opt.lon)};

%% ---- coverage report ----------------------------------------------
info.cover = struct('T', mean(isfinite(F.T(:))), 'S', mean(isfinite(F.S(:))), ...
                    'u', mean(isfinite(F.u(:))), 'N2', mean(isfinite(F.N2(:))));
% Where each variable's usable range actually starts and stops on this grid.
% For a mode fit this matters more than the overall coverage fraction: it is
% the top and bottom edges that bias the fit, not scattered interior gaps.
for f = {'T','S','u','v','N2'}
    good = mean(isfinite(F.(f{1})), 2) > 0.8;      % levels populated >80% of the time
    if any(good)
        info.range.(f{1}) = [z(find(good,1)) z(find(good,1,'last'))];
    else
        info.range.(f{1}) = [NaN NaN];
    end
end
fprintf('  usable depth range (levels >80%% populated in time):\n');
for f = {'T','S','u','v','N2'}
    fprintf('    %-3s %4.0f - %4.0f m\n', f{1}, info.range.(f{1})(1), info.range.(f{1})(2));
end
fprintf(['build_dod_mooring M%d: %d levels x %d hours (%.1f d), %s to %s\n' ...
         '  native %.2f min T/S (k=%d, %d sample dropped by a block mean) | ' ...
         'ADCP %.1f min past the hour, moved by ''%s''\n' ...
         '  coverage  T %.1f%%  S %.1f%%  u %.1f%%  N2 %.1f%%\n'], ...
    m, nz, nh, t_h(end)-t_h(1), datestr(t_h(1)), datestr(t_h(end)), ...
    info.t_step_min, info.avg_k, info.dropped, info.adcp_offset_min, ...
    opt.uv_method, 100*info.cover.T, 100*info.cover.S, 100*info.cover.u, ...
    100*info.cover.N2);
end

%% ======================================================================
function [B, TZ, q] = qc_reject(B, TZ, opt, m)
% Range-reject. Reports per field and per column so a block failure is visible
% as a block rather than as a slightly worse percentage.
lim = struct('Sdbar', [32 35.6], 'Cdbar', [2.5 6.5], ...
             'TSdbar', [3 32],    'PSdbar', [10 1200]);
q = struct('applied', opt.qc, 'dropped', struct(), 'T_dropped', 0);
for f = fieldnames(lim).'                       % so the readme can always read
    q.dropped.(f{1}) = 0;                       % these, qc on or off
end
if ~opt.qc
    q.note = 'opt.qc FALSE: nothing rejected -- M2 C = -25 and the S spikes are still in';
    warning('build_dod_mooring:noqc', '%s', q.note);
    return
end
for f = fieldnames(lim).'
    v = double(B.(f{1}));
    r = lim.(f{1});
    bad = isfinite(v) & (v < r(1) | v > r(2));
    q.dropped.(f{1}) = sum(bad(:));
    if any(bad(:))
        v(bad) = NaN;
        B.(f{1}) = v;
        cols = find(any(bad, 1));
        fprintf('  qc M%d %s: dropped %d (%.4f%%) outside [%g %g], col(s) %s\n', ...
            m, f{1}, sum(bad(:)), 100*mean(bad(:)), r(1), r(2), mat2str(cols));
    end
end
badT = isfinite(TZ) & (TZ < 3 | TZ > 32);       % expected to be a no-op
q.T_dropped = sum(badT(:));
if any(badT(:))
    TZ(badT) = NaN;
    fprintf('  qc M%d TZdepth: dropped %d outside [3 32] degC (unexpected -- look)\n', ...
        m, q.T_dropped);
end
q.note = 'opt.qc true: S [32 35.6], C [2.5 6.5], T [3 32], P [10 1200], before binning';
end

%% ======================================================================
function Y = hbin(t, X, t_h)
% Centred hourly bin mean of X [nrow x nt] onto the exact-hour grid t_h.
%
% Each sample goes to the NEAREST hour, so the bin around t_h(k) is
% [t_h(k)-30min, t_h(k)+30min) -- symmetric, hence no phase shift. A leading
% block mean (samples 1:6, 7:12, ...) is what reshape gives you and it lands
% 25 minutes late, which is 12 deg of M2 phase for free.
%
% Done as one sparse matrix multiply rather than a loop over hours, and it
% counts the finite samples per bin instead of assuming six, so gaps and the
% partial end hour come back NaN rather than biased.
t   = double(t(:)).';
if size(X, 2) ~= numel(t)
    X = X.';                                   % tolerate [nt x nrow]
end
nh  = numel(t_h);
idx = round((t - t_h(1))*24) + 1;              % nearest hour
ok  = idx >= 1 & idx <= nh;
M   = sparse(find(ok), idx(ok), 1, numel(t), nh);

G   = isfinite(X);
X0  = X;  X0(~G) = 0;
num = X0 * M;
den = double(G) * M;
Y   = num ./ den;                              % 0/0 -> NaN, which is right
Y(den == 0) = NaN;
end

%% ======================================================================
function S = s_from_pressure(P_i, S_i, z, extend_deep)
% Interpolate the ns instruments in pressure, hour by hour, honouring the
% time-varying PSdbar.
%
% The two ends are not symmetric. ABOVE the shallowest sensor, S is left NaN:
% the near-surface is where S actually varies (fresh lenses, barrier layers),
% so extrapolating up from 60-94 dbar would invent the one signal the region
% is known for. BELOW the deepest sensor, S is held constant when
% extend_deep -- deep S std at these moorings is 0.004-0.011 psu, so the cost
% is ~0.01 psu, and the alternative is losing 220-320 m of good T out of N2.
nh = size(P_i, 2);  S = nan(numel(z), nh);
for k = 1:nh
    p = P_i(:,k);  s = S_i(:,k);
    g = isfinite(p) & isfinite(s);
    if sum(g) < 2, continue; end
    [p, iu] = unique(p(g));  s = s(g);  s = s(iu);
    S(:,k) = interp1(p, s, z(:), 'linear', NaN);
    if extend_deep
        S(z(:) > p(end), k) = s(end);
    end
end
end

%% ======================================================================
function S = s_from_ts(T_i, S_i, T, mask_gap)
% Build the pooled S(T) curve with ts_relation, then map it through the
% gridded T field. The curve is built by the SAME function plot_dod_TS draws,
% so what is plotted is what is used -- they cannot drift apart.
if nargin < 4 || isempty(mask_gap), mask_gap = true; end
try
    [tc, sc, rel] = ts_relation(T_i, S_i);
catch ME
    warning('build_dod_mooring:tsfit', ...
        'Could not build an S(T) curve (%s); S is NaN.', ME.message);
    S = nan(size(T));  return
end
S = reshape(interp1(tc, sc, T(:), 'linear', NaN), size(T));

% Water inside a void in the T-S data gets no salinity. Across these moorings
% the void is about 6.9-8.0 degC, the 290-720 dbar hole between the 4th and
% deepest instruments. Interpolating across it would read S off a straight
% line drawn between two unrelated water masses and hand N2 a fabricated
% gradient at exactly those depths.
if mask_gap && ~isempty(rel.gap_range)
    for r = 1:size(rel.gap_range, 1)
        lo = rel.gap_range(r,1);  hi = rel.gap_range(r,2);
        S(T > lo & T < hi) = NaN;
        fprintf(['  s_method ts: %.2f-%.2f degC is a void in the T-S data, ' ...
                 'S left NaN there (%d cells)\n'], lo, hi, ...
                 sum(T(:) > lo & T(:) < hi));
    end
end
end

%% ======================================================================
function [N2, note] = n2_of(T, S, z, opt)
% Buoyancy frequency, if gsw is on the path. Returned on the z grid, having
% been computed on the mid-points where the gradient actually lives.
N2 = nan(size(T));
if isempty(which('gsw_Nsquared'))
    note = 'N2: all NaN -- gsw not on the path (addpath /Users/tej/Documents/MATLAB/gsw)';
    warning('build_dod_mooring:nogsw', '%s', note);
    return
end
p   = gsw_p_from_z(-z(:), opt.lat);
SA  = gsw_SA_from_SP(S, repmat(p, 1, size(S,2)), opt.lon, opt.lat);
CT  = gsw_CT_from_t(SA, T, repmat(p, 1, size(T,2)));
[n2m, pmid] = gsw_Nsquared(SA, CT, repmat(p, 1, size(T,2)), opt.lat);
for k = 1:size(T,2)
    g = isfinite(n2m(:,k)) & isfinite(pmid(:,k));
    if sum(g) > 1
        N2(:,k) = interp1(pmid(g,k), n2m(g,k), p, 'linear', NaN);
    end
end
note = 'N2: gsw_Nsquared on mid-points, interpolated back to z; inherits the S approximation';
end

%% ======================================================================
function str = deep_s_note(extend)
if extend
    str = 'held constant (deep S std is 0.004-0.011 psu, so this costs ~0.01 psu)';
else
    str = 'left NaN, which also truncates N2 there and loses 220-320 m of good T';
end
end

%% ======================================================================
function opt = setdef(opt, f, v)
if ~isfield(opt, f) || isempty(opt.(f)), opt.(f) = v; end
end
