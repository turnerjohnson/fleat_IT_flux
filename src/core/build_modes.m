function M = build_modes(zfull, N2full, lat, nbc, rho0)
% BUILD_MODES  Baroclinic displacement, velocity and pressure modes for a flat bottom.
%
%   M = build_modes(zfull, N2full, lat, nbc, rho0)
%
%   zfull   [nz x 1]  depth grid, POSITIVE DOWN, starting at 0, ending at H.
%   N2full  [nz x 1]  buoyancy frequency squared [s^-2], strictly positive, smoothed.
%   lat     central latitude in degrees (used only for local g).
%   nbc     number of BAROCLINIC modes wanted (cfg.nbc, currently 6).
%   rho0    reference density, default 1027 kg m^-3.
%
%   Returns struct M with
%     M.Zeta  [nz x nbc]  displacement modes, baroclinic 1..nbc, max|Zeta| = 1
%     M.Zu    [nz x nbc]  horizontal velocity modes, baroclinic 1..nbc, max|Zu| = 1
%     M.Zbt   [nz x 1]    barotropic velocity mode (essentially a constant)
%     M.Zu0   [nz x nbc+1] the same velocity modes with the BAROTROPIC MODE IN
%                         COLUMN 1, i.e. column j+1 is mode j and column 1 is
%                         mode 0. See INDEXING below.
%     M.P     [nz x nbc]  pressure modes [Pa per metre of displacement amplitude, + up]
%     M.c     [1 x nbc]   eigenspeeds [m/s]
%     M.z, M.N2, M.H, M.rho0
%
%   THREE THINGS ABOUT sw_vmodes THAT WILL BITE YOU IF YOU DO NOT KNOW THEM
%   ----------------------------------------------------------------------
%   1. Its FIRST returned column is the BAROTROPIC mode, not baroclinic mode 1.
%      Verified numerically: PVel(1) = sqrt(g*H) and Hori(:,1) is constant to
%      0.5%. So baroclinic mode j lives in column j+1, and you must ask for
%      nbc+1 modes to get baroclinic 1..nbc. This function does that for you.
%
%   2. It rescales Vert and Hori SEPARATELY so each has max absolute value 1.
%      That throws away the amplitude relationship between them, so you cannot
%      convert a displacement amplitude into a pressure amplitude using Hori
%      alone. This function rebuilds the pressure mode by direct integration
%      instead, which is normalization-proof.
%
%   3. Vert is NOT orthogonal under the natural inner product once the
%      barotropic column is included: cols 1 and 2 correlate at -0.88. The
%      baroclinic columns on their own ARE orthogonal. Displacement has no
%      barotropic mode anyway (Nash et al. 2005 eq. 16 sums from j=1), so this
%      function returns only the baroclinic set for Zeta.
%
%   INDEXING -- WHICH MODE LIVES IN WHICH COLUMN
%   --------------------------------------------
%   Displacement has no barotropic mode (Nash et al. 2005 eq. 16 sums from
%   j = 1) but velocity does (eq. 15 sums from j = 0), so no single column
%   convention can serve both. Rather than leave that offset implicit, this
%   function returns velocity twice under two names:
%
%     M.Zu   column j   = BAROCLINIC mode j.   Pairs 1:1 with M.Zeta and M.P.
%     M.Zu0  column j+1 = mode j, column 1 = mode 0 (barotropic).
%
%   Use M.Zu0 wherever the amplitude vectors are 0-based -- uhat/vhat from
%   fit_uv_modes have mode 0 in ROW 1, so M.Zu0 is the matching design matrix
%   and no +1 is needed. Use M.Zu wherever velocity is paired with pressure,
%   as in the flux cross-integrals, where both sides are baroclinic-only.
%   M.Zbt is kept as the barotropic column on its own.
%
%   SIGN CONVENTION
%   ---------------
%   Everything here is z POSITIVE DOWN, matching sw_vmodes and has
%   appropriate negative signs for eta (which is isopycnal depth 
%   anomaly, positive up). With that choice
%
%       rho'(d,t) = -eta(d,t) * drhobar/dd = -eta * rho0*N2/g
%       p'(d,t)   = rho0 * integral_0^d N2*eta dd_hat,  then demeaned
%
%   The leading minus is folded into M.P, so you can just write
%
%       pprime(:,t) = M.P * ahat(:,t)
%
%   with ahat the displacement amplitudes in metres. Physical check: push the
%   isopycnals down everywhere (eta > 0) and you have removed mass from the
%   column above, so interior p' goes negative. That is what M.P delivers.

if nargin < 5 || isempty(rho0), rho0 = 1027; end

zfull  = zfull(:);
N2full = N2full(:);
H      = zfull(end);

if zfull(1) > 0.01
    error('build_modes:grid','zfull must start at (or essentially at) z = 0.');
end
if any(~isfinite(N2full)) || any(N2full <= 0)
    error('build_modes:N2','N2full must be finite and strictly positive everywhere.');
end

g = 9.78049*(1 + 5.2884e-3*sind(lat)^2 - 5.9e-6*sind(2*lat)^2);  % same g as sw_vmodes

% ask for nbc+1 because column 1 comes back barotropic
[Vert, Hori, ~, PVel] = sw_vmodes(zfull, sqrt(N2full), lat, nbc+1);

% --- confirm column 1 really is the barotropic mode before discarding it -----
bt_speed = sqrt(g*H);
if abs(PVel(1) - bt_speed)/bt_speed > 0.05
    warning('build_modes:notBarotropic', ...
        ['Column 1 of sw_vmodes has c = %.1f m/s but sqrt(gH) = %.1f m/s. ' ...
         'The barotropic assumption may not hold for this N2 -- inspect before trusting.'], ...
         PVel(1), bt_speed);
end

M.Zbt  = Hori(:,1);
M.Zeta = Vert(:,2:nbc+1);
M.Zu   = Hori(:,2:nbc+1);
M.Zu0  = Hori(:,1:nbc+1);            % mode 0 in column 1; see INDEXING above
M.c    = PVel(2:nbc+1);

% --- pressure modes by direct hydrostatic integration ------------------------
% P_j(d) = rho0 * int_0^d N2 * Zeta_j dd, with the depth mean removed so that
% int_0^H P_j dd = 0. That last step is Nash et al. (2005) eq. 4: it is how the
% unknown surface pressure gets pinned down.
M.P = rho0 * cumtrapz(zfull, bsxfun(@times, N2full, M.Zeta));
M.P = bsxfun(@minus, M.P, trapz(zfull, M.P)/H);

% --- cross-check: P_j must be proportional to the velocity mode Zu_j ---------
% Both are, up to a constant, the vertical derivative of the displacement mode.
% If this correlation is not ~1 then N2full and the modes are inconsistent,
% usually because N2full was edited after sw_vmodes was called.
for j = 1:nbc
    r = corr(M.P(:,j), M.Zu(:,j));
    if abs(r) < 0.999
        warning('build_modes:PvsZu', ...
            'Mode %d: corr(P, Zu) = %.4f, expected ~1. Check that N2full matches the modes.', j, r);
    end
end

M.z = zfull;  M.N2 = N2full;  M.H = H;  M.rho0 = rho0;  M.g = g;
end