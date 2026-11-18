function F = compute_tidal_flux(ahat_M2, uhat_M2, vhat_M2, M, cfg, opt)
% COMPUTE_TIDAL_FLUX  Phase-averaged, depth-integrated baroclinic energy flux.
%
%   F = compute_tidal_flux(ahat_M2, uhat_M2, vhat_M2, M, cfg, opt)
%
%   ahat_M2  [np x nt]        DISPLACEMENT amplitudes [m], from fit_eta_modes
%   uhat_M2  [(1+nu) x nt]    velocity amplitudes [m/s], ROW 1 IS BAROTROPIC
%   vhat_M2  same
%   M        from build_modes
%   opt .periods   phase-average window, in M2 periods (default 2)
%       .lowpass   true to low-pass the product instead of boxcar (default false)
%
%   F.Fx, F.Fy       [1 x nt]  phase-averaged flux, W/m (positive = east/north)
%   F.Fx_inst, F.Fy_inst  [1 x nt]  the product BEFORE the phase average. Not
%                    a flux estimate -- see WHY THE AVERAGE IS NOT OPTIONAL.
%   F.mag, F.dir     magnitude [W/m] and direction [deg true, 0 = north]
%   F.Fx_mode        [nflux x nt] per-mode zonal contribution
%   F.C              [nflux x 1]  the cross-integrals, Pa
%   F.nflux          how many modes actually contributed
%
%   ------------------------------------------------------------------
%   THE EQUATION
%   ------------------------------------------------------------------
%   u'(z,t) = sum_i uhat_i(t) Zu_i(z)        Zu dimensionless
%   p'(z,t) = sum_j ahat_j(t) P_j(z)         P in Pa per metre of displacement
%
%   int u'p' dz = sum_i sum_j uhat_i ahat_j int Zu_i P_j dz
%
%   The velocity and pressure modes share vertical structure (both derive from
%   Hori), so int Zu_i P_j dz = 0 for i /= j and the double sum collapses to
%
%       int u'p' dz = sum_i uhat_i(t) ahat_i(t) C_i,   C_i = int Zu_i P_i dz
%
%   C_i is NOT a norm -- Zu and P have different units, so C_i is their cross
%   integral, in Pa. Units then work out as [m/s][m][Pa] = W/m.
%
%   The sum runs to min(np, nu), both read off the arrays you pass in rather
%   than from cfg: a pressure mode with no velocity partner contributes nothing,
%   since you never measured coherent velocity at that wavenumber. This is
%   Nash's low-pass argument made concrete. With one cfg.nbc driving both fits
%   the two counts normally match, but the min is kept so that narrowing the
%   velocity fit (opt.nmodes in fit_uv_modes) stays safe without touching this.
%
%   Both sides here are BAROCLINIC-ONLY, so M.Zu and M.P pair column for
%   column with no offset. The single +1 below is on uhat_M2/vhat_M2, whose
%   ROW 1 is mode 0. See the INDEXING note in build_modes.
%
%   ------------------------------------------------------------------
%   WHY THE AVERAGE IS NOT OPTIONAL
%   ------------------------------------------------------------------
%   For two narrowband M2 signals, their product is
%
%       A cos(wt + phi_u) cos(wt + phi_p)
%           = (A/2)[ cos(phi_u - phi_p) + cos(2wt + phi_u + phi_p) ]
%
%   a constant plus an oscillation at TWICE the tidal frequency. Nash notes the
%   instantaneous u'p' varies as cos^2(wt + phi), so averaging is necessary for
%   the estimate to mean anything. The average kills the 2w term and leaves the
%   flux, which depends on the PHASE DIFFERENCE between velocity and pressure.
%
%   That is the physics of the thing: a propagating wave has u' and p' in phase
%   (cos 0 = 1, maximum flux); a standing wave has them in quadrature
%   (cos 90 = 0, no net flux). Alford's Figures 5 and 6 plot exactly this
%   relative phase for the same reason.
%
%   A boxcar over an integer number of wave periods removes the 2w term
%   exactly. M2 is 12.4206 h and your sampling is hourly, so an integer number
%   of periods is not an integer number of samples -- 2 periods is 24.84 h,
%   rounded to 25. That leaves ~0.6% leakage. More periods reduces it (10
%   periods = 124.2 h, ~0.16%) at the cost of temporal resolution. Setting
%   opt.lowpass uses a filter instead, which sidesteps the rounding entirely.

if nargin < 6, opt = struct(); end
if ~isfield(opt,'periods'), opt.periods = 2;     end
if ~isfield(opt,'lowpass'), opt.lowpass = false; end

nt     = size(ahat_M2, 2);
np     = size(ahat_M2, 1);
nu     = size(uhat_M2, 1) - 1;          % minus the barotropic row
nflux  = min(np, nu);

fprintf('flux: %d pressure modes, %d velocity modes -> summing %d\n', ...
        np, nu, nflux);
if np ~= nu
    fprintf('  (%d mode(s) have a partner on only one side and contribute nothing)\n', ...
            abs(np - nu));
end

%% ---- cross-integrals -------------------------------------------------
F.C = zeros(nflux,1);
for i = 1:nflux
    F.C(i) = trapz(M.z, M.Zu(:,i) .* M.P(:,i));
end

% off-diagonal check: if these are not small, the collapse does not hold and
% the single sum below is missing cross-mode terms
offmax = 0;
for i = 1:nflux
    for j = 1:nflux
        if i ~= j
            offmax = max(offmax, abs(trapz(M.z, M.Zu(:,i).*M.P(:,j))));
        end
    end
end
fprintf('  cross-integrals C_i [Pa]: %s\n', mat2str(round(F.C',2)));
fprintf('  largest off-diagonal / largest diagonal: %.2e\n', offmax/max(abs(F.C)));
if offmax/max(abs(F.C)) > 0.01
    warning('compute_tidal_flux:notdiag', ...
        'Modes are not orthogonal to 1%%. The single-sum collapse is dropping real cross terms.');
end

%% ---- instantaneous flux ---------------------------------------------
F.Fx_mode = zeros(nflux, nt);
F.Fy_mode = zeros(nflux, nt);
for i = 1:nflux
    F.Fx_mode(i,:) = uhat_M2(i+1,:) .* ahat_M2(i,:) * F.C(i);   % i+1 skips barotropic
    F.Fy_mode(i,:) = vhat_M2(i+1,:) .* ahat_M2(i,:) * F.C(i);
end
% NOT 'omitnan'. fit_eta_modes and fit_uv_modes both preallocate nan(nm,nt) and
% assign a timestep only when that timestep's fit succeeds, so a NaN in Fx_mode
% is always a whole COLUMN: it never means "this one mode is missing", it means
% "this hour has no fit at all". Omitting those NaNs sums them to 0 and hands
% the boxcar below a real zero, which is how an hour with no data turns into an
% hour of no flux and biases F toward zero for the ~12 h either side of a gap.
% Leaving the NaN in place lets movmean's own 'omitnan' skip the hour, which is
% what the per-mode averages a few lines down already do.
Fx_inst = sum(F.Fx_mode, 1);
Fy_inst = sum(F.Fy_mode, 1);

% Kept on the struct because F.Fx_mode is overwritten by the phase average a
% few lines down and nothing else in the pipeline can reconstruct this: the
% product needs ahat and uhat, which survive only under opt.keep = 'full'.
% These two are cheap (2 x nt doubles) and they are the only record of what
% the average was taken OVER, which is what a figure needs to show that the
% averaging is doing something rather than hiding something.
%
% NOT A FLUX ESTIMATE. See the header: this is (A/2)[cos(dphi) + cos(2wt+...)],
% so it swings about zero at twice the tidal frequency with an amplitude of
% order the flux itself. Only its average over whole periods means anything.
F.Fx_inst = Fx_inst;
F.Fy_inst = Fy_inst;

%% ---- phase average ---------------------------------------------------
T_M2 = 24/cfg.M2_cpd;                       % 12.4206 h
if opt.lowpass
    % low-pass well below 2*M2 to kill the sum-frequency term
    fs = 24/cfg.dt_hr;
    [b,a] = butter(4, (cfg.M2_cpd)/(fs/2), 'low');
    Fx_inst_f = fillmissing(Fx_inst,'linear','EndValues','nearest');
    Fy_inst_f = fillmissing(Fy_inst,'linear','EndValues','nearest');
    F.Fx = filtfilt(b,a,Fx_inst_f);
    F.Fy = filtfilt(b,a,Fy_inst_f);
    F.Fx(~isfinite(Fx_inst)) = NaN;  F.Fy(~isfinite(Fy_inst)) = NaN;
    fprintf('  phase average: low-pass below %.3f cpd\n', cfg.M2_cpd);
else
    nwin = round(opt.periods * T_M2 / cfg.dt_hr);
    err  = abs(nwin*cfg.dt_hr - opt.periods*T_M2)/(opt.periods*T_M2);
    F.Fx = movmean(Fx_inst, nwin, 'omitnan');
    F.Fy = movmean(Fy_inst, nwin, 'omitnan');
    fprintf('  phase average: %d-h boxcar (%.1f M2 periods, %.2f%% leakage)\n', ...
            nwin, opt.periods, 100*err);
    for i = 1:nflux
        F.Fx_mode(i,:) = movmean(F.Fx_mode(i,:), nwin, 'omitnan');
        F.Fy_mode(i,:) = movmean(F.Fy_mode(i,:), nwin, 'omitnan');
    end
end

%% ---- magnitude and direction -----------------------------------------
F.mag = hypot(F.Fx, F.Fy);
F.dir = mod(90 - atan2d(F.Fy, F.Fx), 360);   % degrees true, 0 = north
F.nflux = nflux;

fprintf('\n  record-mean flux: Fx %+.3f, Fy %+.3f kW/m\n', ...
        mean(F.Fx,'omitnan')/1e3, mean(F.Fy,'omitnan')/1e3);
fprintf('  mean magnitude %.3f kW/m, direction %.0f deg true\n', ...
        mean(F.mag,'omitnan')/1e3, ...
        mod(90 - atan2d(mean(F.Fy,'omitnan'), mean(F.Fx,'omitnan')), 360));
fprintf('  per-mode mean |F| [kW/m]: %s\n', ...
        mat2str(round(mean(hypot(F.Fx_mode,F.Fy_mode),2,'omitnan')'/1e3, 3)));
fprintf('\n  Alford finds O(1 kW/m) at his moorings -- compare.\n');
fprintf('  A record-mean much smaller than the typical magnitude means the\n');
fprintf('  flux reverses direction, i.e. a standing rather than propagating wave.\n\n');
end