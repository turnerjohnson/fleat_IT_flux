function B = f_background_v2(M, cfg)
% F_BACKGROUND_V2  Mean profile and N^2 for a moored T/S array.
%
%   B = f_background_v2(M, cfg)
%
%   M    mooring struct with .depth, .S, .T [nsens x nt], positive down,
%        rows ordered shallow to deep, NaNs allowed. .sigma is used if
%        present and computed from S/T/depth if not.
%   cfg  mooring config. Needs .lat and .H. Everything else has a default,
%        so a config written before this function existed still runs.
%
%   OUT  B.zref     [nz x 1]   reference levels = labels for the isopycnals
%        B.sig_t    [nz x 1]   stabilized mean sigma_0 on zref (target surfaces)
%        B.zmid     [nz-1 x 1] midpoints, where N^2 is actually defined
%        B.N2mid    [nz-1 x 1] N^2 there, rejected but unsmoothed
%        B.N2mids   [nz-1 x 1] after the vertical running mean
%        B.zfull    [nzf x 1]  uniform, 0 to cfg.H inclusive
%        B.sig_full [nzf x 1]  mean sigma_0 on zfull, for plotting
%        B.ztop,.zbot,.covS,.info   coverage and run diagnostics
%
%   TWO N^2 PROFILES, BOTH ALWAYS COMPUTED:
%        B.N2full_sigma [nzf x 1]  from the record-mean S and T, via sw_bfrq
%        B.N2full_grid  [nzf x 1]  from the time median of the supplied F.N2
%                                  field; [] if the struct carries no N2
%        B.zgrid,.N2grid           that median before smoothing, for the plot
%        B.N2full       [nzf x 1]  whichever cfg.N2_source picked. THIS is what
%                                  build_modes reads.
%        B.N2_source               'sigma' or 'gridded'
%
%   The two are genuinely independent estimates and they do not have to agree.
%   Look at the middle panel of PLOT_BACKGROUND, then switch with
%   SET_N2_SOURCE without rebuilding anything.
%
% -------------------------------------------------------------------------
% HOW THIS DIFFERS FROM F2_BACKGROUND, AND WHY
% -------------------------------------------------------------------------
% v1 built sigma_0 on a 2 m grid and differentiated it. That is backwards for
% this array. The density is known at 40 levels spaced 10 to 148 m apart; pchip
% to 2 m adds no information, and GRADIENT then differences it over 4 m. The
% derivative of a shape-preserving interpolant is not smooth -- pchip sets each
% knot slope to a weighted harmonic mean of the neighbouring secants, which is
% biased toward the smaller one, so N^2 dips wherever the level spacing changes.
% On F2 that put kinks at 117, 270, 304, 643, 743, 1389, 1399, 1449 and 1469 m,
% every one of them at a spacing transition.
%
% v2 differences at the levels where the data lives, smooths there, and only
% then interpolates. Four consequences:
%
%   1. N^2 COMES FROM SW_BFRQ, not from d(sigma_0)/dz. sw_bfrq references each
%      pair to its own midpoint pressure (adiabatic leveling), so it is correct
%      at all depths. sigma_0 over a 1500 m column is not: on F2 it reads
%      0.99 of the true N^2 at 100-300 m but 0.84 at 700-1200 m and 0.76 below
%      1200 m. Since c ~ (1/pi) int N dz, that bias lands in the eigenspeeds.
%
%   2. SIG_FULL AND N2FULL ARE NO LONGER DERIVATIVES OF EACH OTHER, on purpose.
%      Integrating the sw_bfrq N^2 back up does not return sigma_0 -- on F2 it
%      drifts +0.17 kg/m^3 by the bottom. That gap is thermobaricity, not an
%      error. So sig_t and sig_full stay sigma_0, because f2_eta tracks surfaces
%      of the sigma_0 that compute_density_mooring put in M.sigma and the labels
%      have to match. N2full is the real stratification. Do not expect
%      gradient(sig_full) to equal N2full; they answer different questions.
%
%      v1's header warned that extrapolating N^2 and integrating can produce a
%      non-monotonic density profile and junk modes. That worry is handled here
%      directly instead: N^2 is rejected where non-positive and floored before
%      it ever reaches build_modes, so the thing build_modes actually requires
%      is enforced rather than inferred.
%
%   3. BAD VALUES ARE REJECTED, NOT CLAMPED. v1's max(N2, N2_floor) pinned 53 of
%      758 grid points -- all of 0-104 m -- at 1e-7, then smoothed log10 across
%      a three-decade cliff, which rings. Here a non-positive or complex N^2
%      becomes NaN and the running mean steps over it using real data on both
%      sides, and nothing is smoothed across the mixed-layer junction at all.
%
%   4. THE SMOOTHER IS IN METRES. v1 used a 9-point loess on a 2 m grid: an 18 m
%      window against levels 100-148 m apart through most of the column, five to
%      eight times too narrow to touch anything. See RUNMEAN_Z below.
%
%
% -------------------------------------------------------------------------
% CHOOSING BETWEEN THE TWO
% -------------------------------------------------------------------------
% cfg.N2_source = 'sigma' is the default because every step of it is visible in
% this file: we pick the levels, we reject, we smooth, and the mixed-layer cap
% is stated rather than inherited. 'gridded' is worth switching to when you want
% an estimate that owes nothing to those choices, or when the sigma branch is
% starved -- a mooring with few salinity-bearing sensors still has F.N2 at every
% instrument.
%
% What they do NOT share: the gridded field reaches above our shallowest
% reference level (87 vs 107 m on F2) using something in place of the missing
% salinity that this pipeline cannot audit, and it carries whatever smoothing
% Simpson's processing already applied. Neither is the truth. If the modes move
% a lot between them, that spread is a result, not a problem to be hidden.
%
% author TEJ, with help from Claude
% October 2026

%% ---- 0. config, with defaults so other moorings drop straight in ----
for r = {'lat','H'}
    if ~isfield(cfg, r{1})
        error('f_background_v2:cfg','cfg.%s is required.', r{1});
    end
end
lat   = cfg.lat;
H     = cfg.H;
mname = pick(cfg, 'mooring',       'mooring');
covmin= pick(cfg, 'sig_cover_min', 0.2);
zmode = lower(pick(cfg, 'zref_mode',  'instrument'));
zmerge= pick(cfg, 'zref_merge',    8);
zdz   = pick(cfg, 'zref_dz',       10);
zfdz  = pick(cfg, 'zfull_dz',      2);
nfloor= pick(cfg, 'N2_floor',      1e-7);
n2ml  = pick(cfg, 'N2_ml',         nfloor);
Lsm   = pick(cfg, 'N2_smooth_m',   50);
nbc   = pick(cfg, 'nbc',           4);

for q = {'depth','S','T'}
    if ~isfield(M, q{1})
        error('f_background_v2:input','M.%s is missing.', q{1});
    end
end
D = M.depth;  S = M.S;  T = M.T;
if isfield(M,'sigma') && ~isempty(M.sigma)
    SG = M.sigma;                                  % keep f2_eta's exact labels
else
    SG = sw_pden(S, T, sw_pres(D, lat), 0) - 1000;
end
mD_all = mean(D, 2, 'omitnan');
for k = 1:size(D,1)-1
    b = isfinite(D(k,:)) & isfinite(D(k+1,:));
    if any(b)
        dzk = median(D(k+1,b) - D(k,b));
    else
        % never co-reporting: the record means are all there is, so tolerate an
        % inversion smaller than the level-merge scale rather than guessing
        dzk = mD_all(k+1) - mD_all(k) + zmerge;
        warning('f_background_v2:order', ...
            ['%s rows %d and %d never report together; their order is taken ' ...
             'from record means with a %.0f m tolerance.'], mname, k, k+1, zmerge);
    end
    if dzk < 0
        error('f_background_v2:order', ...
            ['%s rows %d and %d are inverted: row %d sits %.1f m ABOVE row %d ' ...
             'on the samples where both report. Rows of M.depth must run ' ...
             'shallow to deep.'], mname, k, k+1, k+1, -dzk, k);
    end
end

%% ---- 1. which levels have density, and where ------------------------
% Both S and T, because N^2 needs both
if ~isequal(size(S), size(T), size(D))
    error('f_background_v2:size','M.depth, M.S and M.T must be the same size.');
end
has    = isfinite(S) & isfinite(T) & isfinite(D);
Dm     = D;  Dm(~has) = NaN;
B.ztop = min(Dm, [], 1);
B.zbot = max(Dm, [], 1);
B.covS = mean(has, 2);
mD     = mean(Dm, 2, 'omitnan');
use    = B.covS >= covmin & isfinite(mD);
z      = sort(mD(use));

switch zmode
  case 'instrument', B.zref = z([true; diff(z) > zmerge]);
  case 'uniform',    B.zref = (ceil(z(1)/zdz)*zdz : zdz : floor(z(end)/zdz)*zdz)';
  otherwise, error('f_background_v2:zref_mode','Unknown zref_mode "%s".', zmode);
end
nz = numel(B.zref);

% ---- how many levels nbc modes actually need -------------------------
%
%   fit_eta_modes   ngood < max(minobs, nbc+2)    -> nbc unknowns
%   fit_uv_modes    ngood < max(minobs, nbc+3)    -> nbc+1, mode 0 included
%
% which is the nbc <= n-2 rule in the literature. So nbc+2 is the floor for
% eta: at nz == nbc+2 a fit is possible only in an hour where EVERY level
% reports. We require one more than that, because a record that cannot survive
% a single dropout cannot produce a time series, and warn through the old +4
% band rather than refusing to run.
%
% nz is NOT the same quantity the fit gates test. It is the total level count
% over the record, so it only bounds what any given hour can have; F6 has 463
% hours with fewer than 4 density samples against nz = 7. The honest test of
% whether THESE levels separate THESE modes is cfg.fit_maxcond, which asks it
% properly and returns NaN when the answer is no -- see the "WHY minobs IS NOT
% ENOUGH" note in fit_eta_modes. Counting is the cheap screen, not the check.

assert(nz >= nbc + 3, 'f_background_v2:levels', ...
       ['Only %d reference levels for %d modes. The eta fit needs %d good ' ...
        'levels (nbc+2) in an hour, so %d cannot survive a single dropout. ' ...
        'Lower cfg.nbc for this mooring -- see the nbc note in ' ...
        'fleat_config_base.'], nz, nbc, nbc+2, nz);
if nz < nbc + 5
    warning('f_background_v2:levels_marginal', ...
        ['%s: %d reference levels for %d modes. The eta fit needs %d good ' ...
         'levels in an hour and the velocity fit %d, so hours with any ' ...
         'dropout will return NaN. Check what fraction cfg.fit_maxcond ' ...
         'rejects before quoting a mode-%d amplitude.'], ...
         mname, nz, nbc, nbc+2, nbc+3, nbc);
end

dzr = diff(B.zref);
fprintf('[%s] S and T at %d of %d sensors (cover >= %.2f)\n', ...
        mname, nnz(use), numel(use), covmin);
fprintf('       zref: %d levels, %.0f to %.0f m, spacing %.0f/%.0f/%.0f m (min/median/max)\n', ...
        nz, B.zref(1), B.zref(end), min(dzr), median(dzr), max(dzr));

%% ---- 2. bin-mean S, T and sigma at those levels ---------------------
% v1 binned sigma alone. N^2 needs S and T, so all three are carried.
edges = [B.zref(1) - dzr(1)/2;  (B.zref(1:end-1) + B.zref(2:end))/2;  B.zref(end) + dzr(end)/2];

ok = has & isfinite(SG);            % one mask, so S, T and sigma share samples
b  = discretize(D(ok), edges);
gd = ~isnan(b);
bi = b(gd);
Sa = S(ok);  Ta = T(ok);  Ga = SG(ok);

Sb   = accumarray(bi, Sa(gd), [nz 1], @mean, NaN);
Tb   = accumarray(bi, Ta(gd), [nz 1], @mean, NaN);
sbar = accumarray(bi, Ga(gd), [nz 1], @mean, NaN);
cnt  = accumarray(bi, 1,      [nz 1], @sum,  0);

nfill = nnz(~isfinite(Sb) | ~isfinite(Tb) | ~isfinite(sbar));
Sb    = fillmissing(Sb,   'linear', 'EndValues','nearest');
Tb    = fillmissing(Tb,   'linear', 'EndValues','nearest');
sbar  = fillmissing(sbar, 'linear', 'EndValues','nearest');
if nfill > 0
    fprintf('       %d empty bin(s) filled by interpolation\n', nfill);
end

%% ---- 3. stabilize the sigma_0 labels --------------------------------
% Only the isopycnal labels need to be monotonic. N^2 gets its own test below.
[B.sig_t, sinfo] = stabilizeSigma(sbar, cnt);
fprintf('       stabilize: %d level(s) adjusted, max %.4g kg m^-3\n', ...
        sinfo.nfix, max(sinfo.maxadj));
if max(sinfo.maxadj) > 0.05
    warning('f_background_v2:bigadj', ...
       ['Stabilization moved the RECORD-MEAN profile by %.3g kg/m^3. That is a ' ...
        'sensor depth or calibration problem, not real overturning.'], max(sinfo.maxadj));
end

%% ---- 4. N^2 at the midpoints, where it is actually defined ----------
B.zmid  = (B.zref(1:end-1) + B.zref(2:end))/2;
B.N2mid = sw_bfrq(Sb, Tb, sw_pres(B.zref, lat), lat);

%% ---- 5. reject, do not clamp ----------------------------------------
% Stabilizing sigma_0 does not guarantee a stable profile under the local
% reference pressure, so this still has work to do.
bad = ~isfinite(B.N2mid) | imag(B.N2mid) ~= 0 | real(B.N2mid) <= 0;
B.N2mid = real(B.N2mid);  B.N2mid(bad) = NaN;
if any(bad)
    fprintf('       rejected N^2 <= 0 at %d of %d midpoints: %s m\n', ...
            nnz(bad), nz-1, mat2str(round(B.zmid(bad))'));
end

%% ---- 6. vertical running mean, nsqSolo's pattern, window in metres --
[B.N2mids, isfull] = runmean_z(log10(B.N2mid), B.zmid, Lsm);
B.N2mids = 10.^B.N2mids;
fprintf('       smoothed over %.0f m: %d of %d midpoints had a full window\n', ...
        Lsm, nnz(isfull), nz-1);

still = ~isfinite(B.N2mids);
if any(still)
    B.N2mids = 10.^fillmissing(log10(B.N2mids), 'linear', 'EndValues','nearest');
    fprintf('       %d midpoint(s) still empty after smoothing, interpolated\n', nnz(still));
end

%% ---- 7. onto zfull: mixed-layer cap, measured band, deep hold -------
% Nothing is smoothed across these joins. Each piece asserts one thing and the
% joins stay visible rather than being blended into something that looks
% measured. mld is capped at the shallowest midpoint since above that there is
% no density at all.
% LINSPACE, not 0:zfdz:H. v1's colon grid stops at the last multiple of zfdz
% below H -- with H = 1515 and zfdz = 2 it ends at 1514, and build_modes takes
% H = zfull(end), so the solve quietly ran on a 1514 m column. This lands on H
% exactly and stays uniform; zfdz is honoured to the nearest divisor of H.
B.zfull = linspace(0, H, round(H/zfdz) + 1)';
mld     = min(pick(cfg, 'mld', B.zmid(1)), B.zmid(1));
lg      = nan(size(B.zfull));

ml  = B.zfull <= mld;
bl  = B.zfull > mld & B.zfull < B.zmid(1);
mid = B.zfull >= B.zmid(1) & B.zfull <= B.zmid(end);
dp  = B.zfull > B.zmid(end);

lg(ml)  = log10(n2ml);
if any(bl)      % empty, and interp1 would reject the duplicate X, when mld == zmid(1)
    lg(bl) = interp1([mld; B.zmid(1)], log10([n2ml; B.N2mids(1)]), B.zfull(bl));
end
lg(mid) = interp1(B.zmid, log10(B.N2mids), B.zfull(mid));   % linear, no overshoot
lg(dp)  = log10(B.N2mids(end));                             % hold, assert nothing

B.N2full_sigma = max(10.^lg, nfloor);
fprintf('       mixed layer 0-%.0f m at N^2 = %.2g; below %.0f m held at %.2g\n', ...
        mld, n2ml, B.zmid(end), B.N2mids(end));

%% ---- 7b. the other estimate: the gridded product ---------------------
% F.N2 as supplied, median over the record at the mean sensor depths. Median,
% not mean, because N^2 is noisy and skewed and this is what the sanity panel
% has always plotted.
%
% This path is INDEPENDENT of nearly everything above it. It does not use zref,
% the coverage test, stabilizeSigma or the mixed-layer cap, and it reaches
% sensors we dropped: on F2 it has 43 levels from 87 to 1503 m against our 40
% from 107 to 1496 m, because F.N2 is finite at the 87 m sensor that reports no
% salinity at all. That independence is the whole value of carrying it -- and
% also the reason to be careful, since whatever stands in for salinity up there
% is not something this pipeline can check.
%
% Smoothed with the SAME runmean_z window as the sigma branch so the two are
% compared on equal footing. B.N2grid keeps the unsmoothed median for the plot.
B.N2full_grid = [];  B.zgrid = [];  B.N2grid = [];
if isfield(M,'N2') && ~isempty(M.N2) && isequal(size(M.N2), size(D))
    zg = mean(D, 2, 'omitnan');
    ng = median(M.N2, 2, 'omitnan');
    kg = isfinite(zg) & isfinite(ng) & ng > 0;
    zg = zg(kg);  ng = ng(kg);
    [zg, iu] = unique(zg);  ng = ng(iu);

    if numel(zg) >= 3
        B.zgrid = zg;  B.N2grid = ng;
        ngs = 10.^runmean_z(log10(ng), zg, Lsm);

        lgg   = nan(size(B.zfull));
        above = B.zfull < zg(1);
        band  = B.zfull >= zg(1) & B.zfull <= zg(end);
        below = B.zfull > zg(end);

        % Above the shallowest gridded level there is no estimate at all, so
        % this takes the same neutral cap the sigma branch uses rather than
        % inventing a value. cfg.N2_ml defaults to cfg.N2_floor.
        lgg(above) = log10(n2ml);
        lgg(band)  = interp1(zg, log10(ngs), B.zfull(band));

        % Below the deepest level: a straight line, least squares over the
        % deepest N2_grid_fit_m metres, continued to H and floored.
        % A TWO-POINT slope would be identically zero here. The readme says
        % dTdz is "computed as forward average between successive instruments",
        % so F.N2 repeats its deepest value -- on both F2 and F4 the last two
        % levels are equal to all printed digits. Hence the span.
        fitm = pick(cfg, 'N2_grid_fit_m', 100);
        sel  = zg >= zg(end) - fitm;
        if nnz(sel) >= 2
            pc = polyfit(zg(sel), ngs(sel), 1);
        else
            pc = [0 ngs(end)];
        end
        if any(below)
            lgg(below) = log10(max(polyval(pc, B.zfull(below)), nfloor));
        end

        B.N2full_grid = max(10.^lgg, nfloor);
        fprintf(['       gridded N^2: %d levels, %.0f to %.0f m; above %.0f m set to %.2g, ' ...
                 'below extrapolated at %.2e s^-2/m over %.0f m\n'], ...
                numel(zg), zg(1), zg(end), zg(1), n2ml, pc(1), H - zg(end));
    end
end

%% ---- 7c. which one build_modes gets ---------------------------------
src = lower(pick(cfg, 'N2_source', 'sigma'));
if strcmp(src,'gridded') && isempty(B.N2full_grid)
    warning('f_background_v2:nogrid', ...
        ['cfg.N2_source is "gridded" but no usable N2 field came in on the ' ...
         'mooring struct. Falling back to sigma.']);
    src = 'sigma';
end
B = set_N2_source(B, src);
fprintf('       N2full <- %s\n', B.N2_source);

%% ---- 8. sig_full, for plotting only ---------------------------------
% A sigma_0 profile on the same grid, so the mean-profile panel still works.
% It is NOT the integral of N2full and is not used by build_modes.
sig = nan(size(B.zfull));
in  = B.zfull >= B.zref(1) & B.zfull <= B.zref(end);
sig(in) = interp1(B.zref, B.sig_t, B.zfull(in), 'linear');
sig(B.zfull < B.zref(1))  = B.sig_t(1);
dsdz = (B.sig_t(end) - B.sig_t(end-1)) / dzr(end);
bot  = B.zfull > B.zref(end);
sig(bot) = B.sig_t(end) + dsdz*(B.zfull(bot) - B.zref(end));
B.sig_full = sig;

%% ---- 9. what it cost to average before differencing -----------------
% Shaun differences each profile then averages. We cannot: the modes must be
% time-invariant. But time-averaging at a fixed depth convolves the profile with
% the displacement PDF, which broadens and lowers peaks in N^2. This is the size
% of that smearing, reported rather than hidden.
Pi  = sw_pres(D, lat);
n2i = sw_bfrq(S, T, Pi, lat);
n2i = real(n2i);  n2i(n2i <= 0 | ~isfinite(n2i)) = NaN;
zmi = mean((D(1:end-1,:) + D(2:end,:))/2, 2, 'omitnan');
n2i = mean(n2i, 2, 'omitnan');
gi  = isfinite(zmi) & isfinite(n2i);
[zu, iu] = unique(zmi(gi));  nu = n2i(gi);
if numel(zu) > 2
    ref = interp1(zu, nu(iu), B.zmid, 'linear');
    rat = ref ./ B.N2mids;
    fprintf('       mean-then-difference vs difference-then-mean: median ratio %.2f\n', ...
            median(rat, 'omitnan'));
    B.info.smear_ratio = rat;
end

%% ---- 10. diagnostics -------------------------------------------------
[n2max, imax] = max(B.N2full);
fprintf('       N^2 peak %.3g s^-2 at %.0f m (T_b = %.1f min), deep %.3g s^-2\n', ...
        n2max, B.zfull(imax), 2*pi/sqrt(n2max)/60, B.N2full(end));

B.info.nfix      = sinfo.nfix;
B.info.nreject   = nnz(bad);
B.info.nfull     = nnz(isfull);
B.info.smooth_m  = Lsm;
B.info.mld       = mld;
B.info.bincount  = cnt;
B.info.zref_dz   = dzr;
end


% =========================================================================
function v = pick(cfg, name, default)
%PICK  cfg field if it exists and is non-empty, otherwise the default. Keeps
%      this function runnable against configs written before it existed.
if isfield(cfg, name) && ~isempty(cfg.(name))
    v = cfg.(name);
else
    v = default;
end
end
