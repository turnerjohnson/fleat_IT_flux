function U = fit_uv_modes(F2_h, M, cfg, opt)
% FIT_UV_MODES  u' and v' from moored velocity, with the barotropic tide removed.
%
%   U = fit_uv_modes(F2_h, M, cfg, opt)
%
%   F2_h  struct with .u, .v, .depth, .time   [nsens x nt]
%   M     struct from build_modes -- needs M.Zu0, M.Zu, M.z, M.H
%   cfg   from f2_config
%   opt  .nmodes   baroclinic velocity modes to solve for (default cfg.nbc)
%       .weight    'wkb' (default) | 'none'
%       .minobs    min good levels per timestep (default 8)
%
%   U.uhat, U.vhat        [(1+nm) x nt]  amplitudes. ROW 1 IS MODE 0, BAROTROPIC.
%   U.uhat_M2, U.vhat_M2  same, bandpassed to M2
%   U.up, U.vp            [nz x nt]  baroclinic u', v' on the full mode grid, M2 band
%   U.u_bt, U.v_bt        [1 x nt]   barotropic amplitude, M2 band
%   U.ubar, U.vbar        [nsens x 1] time-mean profile that was removed
%   U.resid_frac          [1 x nt]   residual variance not captured by the fit
%
%   ------------------------------------------------------------------
%   THE MODE SETS DIFFER FROM THE DISPLACEMENT CASE
%   ------------------------------------------------------------------
%   Displacement uses Z^xi, which vanishes at both boundaries, and Nash eq. 16
%   sums from j = 1. Velocity uses Z^u = dZ^xi/dz, and Nash eq. 15 sums from
%   j = 0. Mode 0 is the barotropic mode -- depth-independent flow, the tide
%   sloshing the whole column. Displacement has no barotropic mode; velocity
%   does, and it is usually LARGER than the baroclinic part. So you are
%   extracting a small signal riding on a big one, which is why this is harder
%   than the pressure side.
%
%   ------------------------------------------------------------------
%   WHAT THIS IMPLEMENTS
%   ------------------------------------------------------------------
%   Nash eq. 5:  u'(z,t) = u(z,t) - ubar(z) - u_o(t)
%   Nash eq. 6:  (1/H) int_{-H}^{0} u' dz = 0, which determines u_o
%
%   The three terms are the instantaneous velocity, the time-mean profile
%   (mean flow and mesoscale), and a depth-uniform time-varying part -- the
%   barotropic tide. Here ubar is removed explicitly and u_o is obtained as
%   the mode-0 amplitude of the fit, since solving the least-squares problem
%   with a constant column in the design matrix IS the projection that
%   eq. (6) describes.
%
%   Two routes are computed and compared:
%     A  fit modes 0..nm, then drop mode 0
%     B  remove the depth average over the SAMPLED range, then fit 1..nm
%   These differ because B's depth average is incomplete -- at F2 the array
%   spans 87-1503 m of a 1515 m column, so the top 87 m is missing from the
%   integral that eq. (6) defines over the full depth, and at F6 it spans
%   137-420 m of 434 m. Route A extends to full depth via the modes and is the
%   one to use; B is a consistency check, and U.zB records the span it used.
%
%   This is NOT the fit-vs-filter order question (opt.check_order), and not the
%   p' check that lived in the old scripts as pprime_A vs pprime_B. Three
%   different pairs, three different questions: this one is about the unsampled
%   depth range, the band order is about commuting a bandpass with a
%   least-squares fit, and the p' pair was about M.P against a hydrostatic
%   integral.
%
%   ------------------------------------------------------------------
%   ASSUMPTIONS, AND WHICH ONES BITE
%   ------------------------------------------------------------------
%   - Flat bottom. The modes assume it (Alford eq. A2). F2 sits in 1515 m near
%     Palau, i.e. on a slope. Alford restricted his moorings to H >= 3200 m
%     partly for this reason. State it.
%   - Same N2 as the displacement fit. Required, or u' and p' are projected
%     onto inconsistent bases and the correlation is meaningless. This means
%     the mixed-layer assumption propagates into the velocity side too.
%   - Isotropic modal structure: u and v share Z^u and differ only in
%     amplitude. Fit each separately against the same design matrix.
%   - Hydrostatic, omega << N. Fine here: M2 at 1.93 cpd against N of ~26 cpd
%     in the deep water and far higher in the thermocline.
%   - Linear. There is an M4 peak at 3.86 cpd in the eta spectra, so
%     nonlinearity is measurable, if small.
%
%   ------------------------------------------------------------------
%   DEGREES OF FREEDOM -- READ THIS BEFORE TRUSTING ANY NUMBER
%   ------------------------------------------------------------------
%   F2_h.u has 43 rows but the velocity was extrapolated from THREE ADCPs
%   (nominal 76, 789, 1496 m; the 88 m unit failed). So the effective number
%   of independent constraints is ~3, not 43, no matter how many rows the
%   array has. Nash notes that formally M instruments resolve M modes and that
%   geometry and spectral shape cut it below that in practice -- with three
%   instruments you can hope for mode 0 plus one or two baroclinic modes, not
%   four.
%
%   This function CANNOT detect the extrapolation; it will happily fit four
%   modes to 43 rows and return tight-looking amplitudes. The residual will be
%   near zero because the input is smooth by construction, not because the fit
%   is good. Identify the real ADCP depths and pass opt.nmodes accordingly.
%
%   The saving grace: your p' table gives 51.6, 37.9, 12.3, 11.1 Pa by mode, so
%   modes 1-2 carry ~79% of the pressure. Since p' acts as a low-pass on the
%   u'p' correlation (Nash: p' is a vertical integral of xi, hence redder by
%   k_z^2), resolving mode 0 cleanly plus modes 1-2 is probably enough for the
%   flux even though displacement wanted four.

if nargin < 4, opt = struct(); end
if ~isfield(opt,'nmodes'), opt.nmodes = cfg.nbc; end
if ~isfield(opt,'weight'), opt.weight = cfg.fit_weight; end
if ~isfield(opt,'minobs'), opt.minobs = 8; end

nm = opt.nmodes;
if nm > size(M.Zu,2)
    error('fit_uv_modes:nmodes','Asked for %d baroclinic modes, M.Zu has %d.', nm, size(M.Zu,2));
end
if ~isfield(M,'Zu0')
    error('fit_uv_modes:stalemodes', ...
        ['M has no Zu0 field, so it was built before the single-nbc convention. ' ...
         'Re-run build_modes (cell 4) -- a modes struct left in the workspace ' ...
         'from an earlier run is stale.']);
end
if size(M.Zu0,2) < nm+1
    error('fit_uv_modes:nmodes0','M.Zu0 has %d columns, need %d (mode 0 plus %d baroclinic).', ...
          size(M.Zu0,2), nm+1, nm);
end

u  = F2_h.u;   v = F2_h.v;   D = F2_h.depth;   t = F2_h.time(:).';
[ns, nt] = size(u);

fprintf('velocity: %d rows, %d timesteps, %.1f%% u finite\n', ns, nt, 100*mean(isfinite(u(:))));

%% ---- 1. remove the time-mean profile (Nash eq. 5, second term) -------
% Mean flow and mesoscale. The M2 bandpass would remove these anyway, but
% taking them out first improves the conditioning of the fit.
U.ubar = mean(u, 2, 'omitnan');
U.vbar = mean(v, 2, 'omitnan');
up0 = u - U.ubar;
vp0 = v - U.vbar;

%% ---- 2. WKB weights, as for the displacement fit ---------------------
zmean = mean(D, 2, 'omitnan');
switch lower(opt.weight)
  case 'none'
    w = ones(ns,1);
  case 'wkb'
    Nz   = sqrt(M.N2);
    Nbar = trapz(M.z, Nz)/M.H;
    zp   = cumtrapz(M.z, Nz/Nbar);
    zpr  = interp1(M.z, zp, min(max(zmean, M.z(1)), M.z(end)));
    [zs, o] = sort(zpr);
    edg = [zs(1); (zs(1:end-1)+zs(2:end))/2; zs(end)];
    wsorted = diff(edg);
    w = zeros(ns,1); w(o) = wsorted;
    w = w/mean(w(w>0));
    w(w<=0) = eps;
  otherwise
    error('fit_uv_modes:weight','Unknown weight "%s".', opt.weight);
end

%% ---- 3. fit, timestep by timestep -----------------------------------
% The design matrix is rebuilt each timestep because the sensors MOVE -- 56 m
% of median knockdown at F2. Evaluating the modes at each sensor's
% instantaneous depth is the correct thing; using mean depths would smear the
% mode structure by the knockdown amplitude.
Zall = M.Zu0(:,1:nm+1);                  % [nz x (1+nm)], column 1 = mode 0 (barotropic)

U.uhat = nan(1+nm, nt);
U.vhat = nan(1+nm, nt);
U.resid_frac = nan(1, nt);
U.cond = nan(1, nt);
U.w    = w;              % the sweep refits on exactly these weights
if ~isfield(opt,'maxcond') || isempty(opt.maxcond), opt.maxcond = 100; end
nfew = 0;  nill = 0;   % see the note on conditioning in fit_eta_modes

for it = 1:nt
    zi = D(:,it);
    ok = isfinite(zi) & isfinite(up0(:,it)) & isfinite(vp0(:,it));
    if nnz(ok) < max(opt.minobs, nm+3), nfew = nfew + 1; continue; end

    A = interp1(M.z, Zall, min(max(zi(ok), M.z(1)), M.z(end))); % design matrix: modes interp to sensor depths
    if any(~isfinite(A(:))), continue; end

    Wk = sqrt(w(ok));
    Aw = bsxfun(@times, A, Wk);

    % A count of sensors says nothing about whether they span the column.
    % Same failure and same fix as the displacement fit -- see fit_eta_modes.
    U.cond(it) = cond(Aw);
    if U.cond(it) > opt.maxcond, nill = nill + 1; continue; end

    au = Aw \ (up0(ok,it).*Wk); % amplitudes of u
    av = Aw \ (vp0(ok,it).*Wk); % amplitudes of v
    U.uhat(:,it) = au;
    U.vhat(:,it) = av;

    r  = [up0(ok,it); vp0(ok,it)] - [A*au; A*av];
    s  = [up0(ok,it); vp0(ok,it)];
    U.resid_frac(it) = sum(r.^2)/max(sum(s.^2), eps);
end

nfit = nnz(all(isfinite(U.uhat),1));
fprintf('complete fits at %d of %d timesteps (%d too few sensors, %d ill-conditioned)\n', ...
        nfit, nt, nfew, nill);
fprintf('median residual fraction: %.3f\n', median(U.resid_frac,'omitnan'));
fprintf('  (a near-zero residual on extrapolated input means the input was\n');
fprintf('   smooth, not that the fit resolved anything)\n');

%% ---- 4. bandpass the amplitudes -------------------------------------
gap = ~all(isfinite(U.uhat),1);
uf  = fillmissing(U.uhat, 'linear', 2, 'EndValues','nearest');
vf  = fillmissing(U.vhat, 'linear', 2, 'EndValues','nearest');

U.uhat_M2 = bandpass_alford(uf, cfg.dt_hr, cfg.M2_cpd, cfg.bw_c);
U.vhat_M2 = bandpass_alford(vf, cfg.dt_hr, cfg.M2_cpd, cfg.bw_c);
U.uhat_M2(:,gap) = NaN;  U.vhat_M2(:,gap) = NaN;

U.u_bt = U.uhat_M2(1,:);
U.v_bt = U.vhat_M2(1,:);

%% ---- 5. baroclinic u', v' on the full mode grid ---------------------
% Drop mode 0. What remains satisfies (1/H) int u' dz = 0 by construction,
% because the baroclinic modes are orthogonal to the constant over full depth.
U.up = M.Zu(:,1:nm) * U.uhat_M2(2:end,:);
U.vp = M.Zu(:,1:nm) * U.vhat_M2(2:end,:);

U.bt_check = max(abs(trapz(M.z, U.up)/M.H));
fprintf('\ndepth-mean of u'' (Nash eq. 6, want 0): max %.3g m/s\n', U.bt_check);

fprintf('\nrms M2 amplitude (m/s):\n');
fprintf('  mode 0 (barotropic): u %.4f, v %.4f\n', ...
        std(U.uhat_M2(1,:),0,2,'omitnan'), std(U.vhat_M2(1,:),0,2,'omitnan'));
for j = 1:nm
    fprintf('  mode %d:              u %.4f, v %.4f\n', j, ...
            std(U.uhat_M2(j+1,:),0,2,'omitnan'), std(U.vhat_M2(j+1,:),0,2,'omitnan'));
end

bt = std(U.uhat_M2(1,:),0,2,'omitnan');
bc = sqrt(sum(std(U.uhat_M2(2:end,:),0,2,'omitnan').^2));
fprintf('\nbarotropic / baroclinic rms ratio (u): %.2f\n', bt/max(bc,eps));
fprintf('  A ratio well above 1 is normal and is exactly why mode 0 must be\n');
fprintf('  in the fit -- leaving it out puts that energy into modes 1-2.\n\n');

%% ---- 6. route B, as a cross-check -----------------------------------
% Remove the depth average over the SAMPLED range only, then fit baroclinic
% modes. This is incomplete, since the array never spans the full column -- at
% F2 it misses the top 87 m of 1515 m, at F6 the top 137 m and the bottom 14 m
% of 434 m -- so disagreement with route A localises how much the missing
% surface layer matters for the barotropic estimate. The printout below
% measures the span rather than quoting one mooring's.
uB = nan(1+0+nm, nt);
zB = nan(2, nt);                   % the range route B actually integrated over
for it = 1:nt
    zi = D(:,it);
    ok = isfinite(zi) & isfinite(up0(:,it));
    if nnz(ok) < max(opt.minobs, nm+3), continue; end
    [zs, o] = sort(zi(ok));  us = up0(ok,it); us = us(o);
    ubar_t = trapz(zs, us)/(zs(end)-zs(1));       % over the sampled range only
    A = interp1(M.z, M.Zu(:,1:nm), min(max(zs, M.z(1)), M.z(end)));
    Wk = sqrt(w(ok)); Wk = Wk(o);
    uB(2:end,it) = bsxfun(@times,A,Wk) \ ((us - ubar_t).*Wk);
    uB(1,it) = ubar_t;
    zB(:,it) = [zs(1); zs(end)];
end
U.uhat_B = uB;
U.zB     = zB;

% The span is MEASURED, not quoted. This message read "only 87-1503 m" until
% 30 Oct 2026 -- F2's velocity range, hardcoded, with only the column depth
% taken from M.H. Every other mooring printed a sentence about a depth range it
% does not have: F6 samples 137-420 m of 434 m and was told it sampled 87-1503.
% Knockdown moves the range hour to hour, so report the median.
ztop = median(zB(1,:), 'omitnan');
zbot = median(zB(2,:), 'omitnan');
d = std(U.uhat(1,:) - uB(1,:), 0, 2, 'omitnan');
fprintf('barotropic estimate, route A vs B: rms difference %.4f m/s\n', d);
if isfinite(ztop)
    miss = 100*(ztop + (M.H - zbot))/M.H;
    fprintf(['  (B integrates only %.0f-%.0f m of a %d m column -- %.0f%% of ' ...
             'the depth\n'], ztop, zbot, round(M.H), miss);
    fprintf(['   is outside it -- so this gap measures how much the unsampled\n' ...
             '   surface layer matters. Compare it with the mode-0 rms above.)\n\n']);
else
    fprintf('  (no timestep had enough good levels for route B)\n\n');
end
end