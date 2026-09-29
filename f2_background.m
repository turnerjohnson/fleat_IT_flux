function B = f2_background(Dh, SG, cfg)
% F2_BACKGROUND  Mean density profile, mixed-layer cap, and N2 on a full-depth grid.
%
%   B = f2_background(Dh, SG, cfg)
%
%   Dh  [nsens x nt]  instrument depths, positive down, NaNs allowed
%   SG  [nsens x nt]  sigma-theta at those depths
%   cfg               from f2_config
%
%   B.zref    [nz x 1]  reference levels (also the labels for the isopycnals)
%   B.sig_t   [nz x 1]  stabilized mean sigma-theta on zref -- the target surfaces
%   B.zfull   [nzf x 1] 0 : cfg.zfull_dz : H
%   B.sig_full[nzf x 1] mean profile on the full grid, mixed layer included
%   B.N2full  [nzf x 1] smoothed, strictly positive, consistent with sig_full
%   B.ztop, B.zbot      shallowest/deepest sensor at each time
%
%   WHAT THIS IMPLEMENTS
%   Nash et al. (2005) eq. (2): rho'(z,t) = rho(z,t) - rhobar(z), where rhobar is
%   "vertical profiles of an undisturbed water column in the absence of internal
%   waves". Nash notes the background evolves due to mesoscale processes and that
%   for coarsely sampled data it is hard to separate the background state from a
%   wave-perturbed one. This function builds rhobar once, from the record mean.
%   The modes must be time-invariant for the decomposition to mean anything, so
%   this is the right place for a single record-long mean -- the 25-day running
%   mean belongs downstream, in the displacement calculation, not here.
%
%   THE MIXED LAYER
%   The shallowest T sensor is at 78 m. Everything above it is asserted. Alford's
%   recommendation is to declare a mixed layer there, and this function does that
%   in DENSITY space rather than in N2 space. That ordering matters: if you
%   extrapolate N2 directly and then integrate, nothing guarantees the implied
%   density profile is monotonic, and sw_vmodes will quietly produce garbage
%   modes. Building sigma first and differentiating it guarantees consistency.
%
%   Physically, a neutrally stratified cap forces the displacement modes to be
%   linear in z through the mixed layer (since zeta_zz = -(N2/c2) zeta -> 0
%   there), rising from zero at the surface. That is the correct behaviour and it
%   is why zeta does not simply vanish over the top 78 m.

g = 9.78049*(1 + 5.2884e-3*sind(cfg.lat)^2 - 5.9e-6*sind(2*cfg.lat)^2);  % same g as sw_vmodes


%% ---- 1-2. reference levels where DENSITY is measured ---------------
has = isfinite(SG) & isfinite(Dh);
Dm  = Dh;  Dm(~has) = NaN;
B.ztop = min(Dm, [], 1);                 % shallowest density sensor, per hour
B.zbot = max(Dm, [], 1);
B.covS = mean(has, 2);                    % fraction of hours with sigma, per sensor
mD  = mean(Dm, 2, 'omitnan');
use = B.covS >= cfg.sig_cover_min & isfinite(mD);
z   = sort(mD(use));
switch lower(cfg.zref_mode)
  case 'instrument'
    B.zref = z([true; diff(z) > cfg.zref_merge]);
  case 'uniform'
    B.zref = (ceil(z(1)/cfg.zref_dz)*cfg.zref_dz : cfg.zref_dz : ...
              floor(z(end)/cfg.zref_dz)*cfg.zref_dz)';
  otherwise
    error('f2_background:zref_mode','Unknown zref_mode "%s".', cfg.zref_mode);
end
fprintf('density at %d of %d sensors (cover >= %.2f); top level %.0f m\n', ...
        nnz(use), numel(use), cfg.sig_cover_min, B.zref(1));
nz = numel(B.zref);
assert(nz > cfg.nbc + 4, 'Only %d reference levels for %d modes.', nz, cfg.nbc);
fprintf('zref: %d levels, %.0f to %.0f m (mode "%s")\n', ...
        nz, B.zref(1), B.zref(end), cfg.zref_mode);
%% ---- 3. binned mean profile ----------------------------------------
edges = [B.zref(1) - (B.zref(2)-B.zref(1))/2;
         (B.zref(1:end-1) + B.zref(2:end))/2;
         B.zref(end) + (B.zref(end)-B.zref(end-1))/2];

za = Dh(:); sa = SG(:);
ok = isfinite(za) & isfinite(sa);
za = za(ok);  sa = sa(ok);              % MATLAB will not chain sa(ok)(gd)
b  = discretize(za, edges);
gd = ~isnan(b);
sbar = accumarray(b(gd), sa(gd), [nz 1], @mean, NaN);
cnt  = accumarray(b(gd), 1,      [nz 1], @sum,  0);
sbar = fillmissing(sbar, 'linear', 'EndValues','nearest');

%% ---- 4. stabilize ---------------------------------------------------
% Weighted by bin count so poorly sampled levels absorb more of the adjustment.
[B.sig_t, info] = stabilizeSigma(sbar, cnt);
fprintf('stabilize: %d level(s) adjusted, max %.4g, rms %.4g kg m^-3\n', ...
        info.nfix, max(info.maxadj), max(info.rmsadj));
if max(info.maxadj) > 0.05
    warning('f2_background:bigadj', ...
        ['Stabilization moved the mean profile by %.3g kg/m^3. That is large for a ' ...
         'RECORD MEAN and suggests a sensor depth or calibration problem, not ' ...
         'real overturning.'], max(info.maxadj));
end

%% ---- 5. extend to full depth, with a mixed layer on top ------------
B.zfull = (0 : cfg.zfull_dz : cfg.H)';
sig = nan(size(B.zfull));

in  = B.zfull >= B.zref(1) & B.zfull <= B.zref(end);
sig(in) = interp1(B.zref, B.sig_t, B.zfull(in), 'pchip');

% --- above the top sensor: mixed layer, then blend ---
mld = min(cfg.mld, B.zref(1));
top = B.zfull < B.zref(1);
if mld >= B.zref(1) - 1e-6
    % Mixed layer reaches the top sensor: the whole unmeasured cap is neutral.
    sig(top) = B.sig_t(1);
    fprintf('mixed layer: neutral from 0 to %.0f m (the top zref level)\n', B.zref(1));
else
    % Shallower mixed layer than the top sensor. Neutral to mld, then a
    % stratified blend up to the first measured value. sig_ml is set below
    % sig_t(1) by continuing the shallowest observed gradient across the gap,
    % so the blend layer carries realistic stratification instead of a step.
    ml = B.zfull <= mld;
    bl = top & ~ml;
    dsdz_top = (B.sig_t(2)-B.sig_t(1))/(B.zref(2)-B.zref(1));
    sig_ml   = B.sig_t(1) - dsdz_top*(B.zref(1)-mld);
    sig(ml)  = sig_ml;
    sig(bl)  = interp1([mld; B.zref(1)], [sig_ml; B.sig_t(1)], B.zfull(bl));
    fprintf('mixed layer: neutral 0-%.0f m at sigma %.3f, blended to %.3f at %.0f m\n', ...
            mld, sig_ml, B.sig_t(1), B.zref(1));
end

% --- below the bottom sensor: continue the deepest observed gradient ---
bot = B.zfull > B.zref(end);
if any(bot)
    dsdz_deep = (B.sig_t(end)-B.sig_t(end-1))/(B.zref(end)-B.zref(end-1));
    sig(bot)  = B.sig_t(end) + dsdz_deep*(B.zfull(bot)-B.zref(end));
    fprintf('below %.0f m: extended %.0f m at the deepest observed gradient\n', ...
            B.zref(end), cfg.H - B.zref(end));
end
B.sig_full = sig;

%% ---- 6. N2 ----------------------------------------------------------
N2 = (g/cfg.rho0) * gradient(B.sig_full, B.zfull);

% smooth in log space; clamp first so log is defined
N2 = max(N2, cfg.N2_floor);
N2 = 10.^smoothdata(log10(N2), 'loess', cfg.N2_smooth_span);
B.N2full = max(N2, cfg.N2_floor);

[n2max, imax] = max(B.N2full);
fprintf('N2: peak %.3g s^-2 at %.0f m, deep value %.3g s^-2\n', ...
        n2max, B.zfull(imax), B.N2full(end));
fprintf('    buoyancy period at peak: %.1f min\n', 2*pi/sqrt(n2max)/60);
end