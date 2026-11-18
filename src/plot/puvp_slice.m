function S = puvp_slice(out, when, varargin)
% PUVP_SLICE  One timestep of a mooring: u', v', p', u'p' profiles and the flux.
%
%   S = puvp_slice(out, when)
%   S = puvp_slice(out, when, 'obs', true, 'nmode_show', 3, 'figure', 21)
%   S = puvp_slice(out, when, 'periods', 10)     % wider phase average
%
%   out    the struct from fleat_IT_flux RUN WITH opt.keep = 'full'. Everything
%          this needs is already on it -- see WHERE EACH PANEL COMES FROM.
%   when   a datenum, or an integer index into out.time. The nearest hour is
%          used and the one actually plotted comes back as S.it / S.time.
%
%   Five panels, left to right:
%     1  u'(z)      baroclinic zonal velocity, M2 band
%     2  v'(z)      baroclinic meridional velocity, M2 band
%     3  p'(z)      perturbation pressure
%     4  u'p'(z)    depth-resolved flux density, this hour and phase-averaged
%     5  F          the depth-integrated flux vector at this timestep
%
%   Panel 4 is the bridge between 1-3 and 5: integrate the heavy curves over
%   depth and you have the black arrow. eta is still computed and returned on S,
%   it just has no panel -- it carries the same amplitudes as p' and so says
%   nothing about the flux that p' does not.
%
%   ------------------------------------------------------------------
%   WHERE EACH PANEL COMES FROM  (nothing here is recomputed from raw data)
%   ------------------------------------------------------------------
%     u'(z,t) = out.uv.up(:,it)              = M.Zu * uhat_M2(2:end,it)
%     v'(z,t) = out.uv.vp(:,it)              row 1 of uhat_M2 is BAROTROPIC and
%                                            is excluded, by construction
%     eta(z,t)= M.Zeta * out.ahat.band(:,it)
%     p'(z,t) = out.pprime(:,it)             = M.P * out.ahat.band(:,it)
%     F(t)    = out.flux.Fx(it), out.flux.Fy(it)
%
%   Every profile lives on the full mode grid M.z, so they are the fitted field
%   extended to full depth -- including the top 87 m and the bottom 12 m that no
%   sensor ever saw. The circles are the observations the fit actually had, at
%   their instantaneous (knocked-down) depths. Where the line has no circles
%   near it, it is modal extrapolation: at F2 that is the top of panels 1 and 2,
%   and it propagates straight into panel 4.
%
%   ------------------------------------------------------------------
%   THE OBSERVED OVERLAY IS NOT THE PIPELINE'S OWN ORDER OF OPERATIONS
%   ------------------------------------------------------------------
%   The pipeline always fits first and bandpasses the modal amplitudes; see the
%   top of cell 5 in fleat_IT_flux for why. To draw observations in the same
%   band this has to bandpass
%   the SENSOR series instead, which is the other order. fleat_IT_flux's
%   opt.check_order measures what that costs -- 8-14% at F2 -- so expect the
%   circles to sit near the line, not exactly on it. A systematic offset is the
%   order of operations; a wrong SHAPE is the fit failing.
%
%   For velocity the overlay additionally has to take the barotropic tide back
%   out of the observation, since u'(z) is baroclinic-only:
%       u'_obs(z,t) = bandpass(u - ubar) - uhat_M2(1,t)*Zbt(z)
%   Pass 'obs', false to skip the overlay and the two bandpasses it costs.
%
%   ------------------------------------------------------------------
%   PANEL 4: WHY THE FAINT CURVES DO NOT CARRY THE SIGN
%   ------------------------------------------------------------------
%   The faint curves are u'(z)p'(z) at this one hour. They are an instantaneous
%   product of two narrowband M2 signals, so at every depth they carry a term at
%   2*omega on top of the mean, and for most of the tidal cycle that term is
%   larger than the mean. The depth integral of a faint curve therefore has
%   whatever sign the hour happened to land on, and stepping one hour forward
%   can flip it. It is not a flux.
%
%   The heavy curves are the same product averaged over the same window
%   compute_tidal_flux used (cfg.flux_periods tidal periods, 25 h by default).
%   That average kills the 2*omega term and leaves <u'p'>(z), which IS a flux
%   density: trapz(S.z, S.upp) returns S.F(1) to round-off, and the check is
%   printed. This is the curve to read for where in the column the flux lives.
%
%   The printed check is therefore a real one, and it tests two things at once:
%   that this window matches the boxcar compute_tidal_flux used, and that the
%   modes are orthogonal enough for its single-sum collapse to hold (S.F_depth
%   integrates the reconstructed profiles and so carries the cross-mode terms
%   that S.F drops). Either failing shows up as a percentage.
%
%   ------------------------------------------------------------------
%   PANEL 5: TWO DIFFERENT THINGS CALLED F
%   ------------------------------------------------------------------
%   The heavy black arrow is out.flux, which is phase-averaged over
%   cfg.flux_periods tidal periods (25 h by default). That is the flux.
%
%   The grey dashed arrow is the UNAVERAGED depth integral at this one hour,
%   sum_i uhat_i(t) ahat_i(t) C_i. It is not an estimate of anything on its own:
%   the product of two narrowband M2 signals carries a term at 2*omega, so this
%   arrow swings through a full cycle twice per tide while the black one barely
%   moves. Watching the two separate as you step through hours is the clearest
%   demonstration of why compute_tidal_flux averages at all.
%
%   The thin coloured arrows are the per-mode phase-averaged contributions,
%   out.flux.Fx_mode/Fy_mode, which sum to the black arrow.
%
%   S  .it .time                  the timestep actually drawn
%      .z                         full mode grid [nz x 1]
%      .up .vp .pp .eta           the profiles [nz x 1]; eta is not plotted
%      .pp_mode                   [nz x nflux] per-mode pressure contributions
%      .upp_inst .vpp_inst        u'p'(z) at this hour alone [W/m^2]
%      .upp .vpp                  <u'p'>(z), averaged over S.iwin [W/m^2]
%      .iwin                      the hours in the phase-average window
%      .zobs .uobs .vobs          observed velocity overlay (NaN if 'obs' false)
%      .zref .etaobs              observed eta overlay
%      .F .F_inst .F_mean         [Fx Fy] in W/m, phase-averaged / raw / record
%      .F_mode                    [nflux x 2]
%      .F_depth                   [Fx Fy] from trapz(z, <u'p'>) -- panel 4
%                                 integrated, should equal .F
%      .F_trapz                   [Fx Fy] from trapz(z, u'.*p') at this hour --
%                                 the same number as .F_inst by another route
%      .ax                        the five axes
%
%   S.F_trapz integrates the reconstructed profiles directly instead of using
%   the modal cross-integrals C_i. The two agree to round-off when the modes are
%   orthogonal, which is the assumption compute_tidal_flux's single-sum collapse
%   rests on, so the difference is printed as a check.
%
%   author TEJ, with help from Claude
%   20 October 2026

%% ==== 0. arguments and the guard that matters ========================
ip = inputParser;
ip.addParameter('obs',        true);    % overlay the observations
ip.addParameter('nmode_show', 3);       % per-mode arrows to draw in panel 5
ip.addParameter('periods',    []);      % phase-average window, tidal periods
ip.addParameter('figure',     []);
ip.addParameter('verbose',    true);
ip.parse(varargin{:});
o = ip.Results;

need = {'Mh','modes','ahat','uv','pprime','flux','time'};
miss = need(~isfield(out, need));
if ~isempty(miss)
    error('puvp_slice:notfull', ...
        ['out is missing %s, so it was built with opt.keep = ''flux''.\n' ...
         'Re-run with opt.keep = ''full'' -- the profiles are not stored\n' ...
         'otherwise and cannot be rebuilt from the flux alone.'], ...
         strjoin(miss, ', '));
end

M    = out.modes;
cfg  = out.cfg;
ahat = out.ahat.band;
UV   = out.uv;
t    = out.time(:).';
nt   = numel(t);

% The window panel 4 averages over has to be the one compute_tidal_flux used or
% the heavy curves will not integrate to the arrow. fleat_IT_flux writes the
% value it resolved back onto cfg, so this reads what actually ran rather than
% what the config file says.
if isempty(o.periods)
    if isfield(cfg,'flux_periods') && ~isempty(cfg.flux_periods)
        o.periods = cfg.flux_periods;
    else
        o.periods = 2;
    end
end
if isfield(cfg,'flux_lowpass') && isequal(cfg.flux_lowpass, true)
    warning('puvp_slice:lowpass', ...
        ['This run phase-averaged with opt.flux_lowpass, a Butterworth filter ' ...
         'with no finite window. Panel 4 falls back to a %g-period boxcar, so ' ...
         'its depth integral will not land exactly on out.flux.'], o.periods);
end

%% ==== 1. which timestep ==============================================
if ~isscalar(when) || ~isfinite(when)
    error('puvp_slice:when', 'when must be a scalar datenum or index.');
end
if when > 7e5                              % a datenum, not an index
    [dt_, it] = min(abs(t - when));
    if dt_ > 1                             % more than a day off the record
        warning('puvp_slice:faroff', ...
            'Nearest hour to %s is %s, %.1f days away.', ...
            datestr(when), datestr(t(it)), dt_);
    end
else
    it = round(when);
    if it < 1 || it > nt
        error('puvp_slice:index', 'Index %d is outside 1..%d.', it, nt);
    end
end
S.it = it;  S.time = t(it);

if ~all(isfinite(ahat(:,it)))
    warning('puvp_slice:nofit', ...
        ['No complete modal fit at %s (index %d): the displacement amplitudes ' ...
         'are NaN there, so the profiles will be empty. Pick an hour inside ' ...
         'the complete-fit window.'], datestr(t(it)), it);
end

%% ==== 2. the profiles =================================================
% Straight off out. The only arithmetic is eta, which the pipeline never needs
% on the full grid and so does not store; it has no panel but is cheap and the
% caller may want it.
S.z   = M.z(:);
S.up  = UV.up(:,it);
S.vp  = UV.vp(:,it);
S.eta = M.Zeta * ahat(:,it);
S.pp  = out.pprime(:,it);

nflux   = out.flux.nflux;
nshow   = min(o.nmode_show, nflux);
S.pp_mode = M.P(:,1:nflux) .* ahat(1:nflux,it).';     % [nz x nflux]

%% ==== 3. the observations the fit actually saw =======================
Mh = out.Mh;
S.zref   = Mh.bkgd.zref(:);
S.etaobs = nan(size(S.zref));
S.zobs   = nan(size(Mh.depth,1), 1);
S.uobs   = S.zobs;  S.vobs = S.zobs;

if o.obs
    % eta in the band, at the reference levels. The pipeline keeps this on
    % out.eta.band, which is the series the fit itself was built from under
    % 'filter_first'; re-derive it only for an out saved before that field
    % existed.
    if isfield(out.eta, 'band')
        etab = out.eta.band;
    else
        etab = band_rows(out.eta.eta, cfg);
    end
    S.etaobs = etab(:,it);

    % velocity: mean profile out (as fit_uv_modes does), band, then the
    % barotropic mode out, so what is left is comparable to up/vp.
    ub = band_rows(Mh.u - UV.ubar, cfg);
    vb = band_rows(Mh.v - UV.vbar, cfg);
    zi = Mh.depth(:,it);
    zc = min(max(zi, M.z(1)), M.z(end));            % clamp, as the fits do
    Zbt = interp1(M.z, M.Zbt, zc);
    S.zobs = zi;
    S.uobs = ub(:,it) - UV.uhat_M2(1,it)*Zbt;
    S.vobs = vb(:,it) - UV.vhat_M2(1,it)*Zbt;
end

%% ==== 4. the flux at this hour =======================================
C = out.flux.C(:);
S.F      = [out.flux.Fx(it)  out.flux.Fy(it)];                   % phase-averaged
S.F_mode = [out.flux.Fx_mode(:,it)  out.flux.Fy_mode(:,it)];     % per mode
S.F_mean = [mean(out.flux.Fx,'omitnan')  mean(out.flux.Fy,'omitnan')];

% the raw depth integral at this one hour, no averaging
ui = UV.uhat_M2(2:nflux+1, it);
vi = UV.vhat_M2(2:nflux+1, it);
ai = ahat(1:nflux, it);
% No 'omitnan', for the same reason compute_tidal_flux does not use one: a
% NaN here is a whole hour with no fit, not a missing mode, and summing it
% to 0 would draw a zero-length 'this hour' arrow where there is no
% estimate at all -- and would disagree with S.F_trapz, which goes NaN.
S.F_inst = [sum(ui.*ai.*C)  sum(vi.*ai.*C)];

% the same thing by direct integration of the profiles drawn in panels 1-3.
% Equal to S.F_inst when the modes are orthogonal; that is the assumption.
S.F_trapz = [trapz(S.z, S.up.*S.pp)  trapz(S.z, S.vp.*S.pp)];

%% ==== 5. depth-resolved flux density, panel 4 ========================
% The faint curves: one hour's product, 2*omega term and all.
S.upp_inst = S.up .* S.pp;                        % [W/m^2]
S.vpp_inst = S.vp .* S.pp;

% The heavy curves: the same product averaged over the window compute_tidal_flux
% used. movmean centres an even window on the current and PREVIOUS elements,
% i.e. [it-k/2, it+k/2-1], and shrinks it at the record ends -- both matched
% here so that trapz(S.z, S.upp) lands on out.flux.Fx(it).
nwin   = round(o.periods * (24/cfg.M2_cpd) / cfg.dt_hr);
iw     = (it - floor(nwin/2)) : (it + ceil(nwin/2) - 1);
S.iwin = iw(iw >= 1 & iw <= nt);

S.upp = mean(UV.up(:,S.iwin) .* out.pprime(:,S.iwin), 2, 'omitnan');
S.vpp = mean(UV.vp(:,S.iwin) .* out.pprime(:,S.iwin), 2, 'omitnan');
S.F_depth = [trapz(S.z, S.upp)  trapz(S.z, S.vpp)];

if o.verbose
    den = max(hypot(S.F_inst(1), S.F_inst(2)), eps);
    fprintf('[%s] %s (index %d of %d)\n', out.name, datestr(S.time), it, nt);
    fprintf('  F phase-averaged : %+8.3f %+8.3f kW/m  |F| %.3f  dir %3.0f deg true\n', ...
            S.F/1e3, hypot(S.F(1),S.F(2))/1e3, out.flux.dir(it));
    fprintf('  F this hour, raw : %+8.3f %+8.3f kW/m  (carries the 2*omega term)\n', ...
            S.F_inst/1e3);
    fprintf('  F record mean    : %+8.3f %+8.3f kW/m\n', S.F_mean/1e3);
    fprintf('  modal sum vs trapz(z, u''p''): %.2e relative\n', ...
            norm(S.F_inst - S.F_trapz)/den);
    % Panel 4 integrated against the arrow it is supposed to produce. A few
    % percent means a gap inside the window; see the note in the header.
    fprintf('  trapz(z, <u''p''>) over %d h: %+8.3f %+8.3f kW/m (%.1f%% off F)\n', ...
            numel(S.iwin), S.F_depth/1e3, ...
            100*norm(S.F_depth - S.F)/max(hypot(S.F(1),S.F(2)), eps));
end

%% ==== 6. draw =========================================================
if isempty(o.figure), fh = figure; else, fh = figure(o.figure); clf(fh); end
set(fh, 'Color','w', 'Position', [60 60 1280 520]);
cols = orderedcolors("gem");
zl   = [0 M.H];

ax = gobjects(1,5);

% ---- 1-2: u', v' -----------------------------------------------------
specs = { S.up, S.uobs, 'u'' (m s^{-1})', 'zonal velocity'
          S.vp, S.vobs, 'v'' (m s^{-1})', 'meridional velocity' };
for k = 1:2
    ax(k) = subplot(1,5,k);
    plot(specs{k,1}, S.z, 'k-', 'LineWidth', 1.8); hold on
    if o.obs
        plot(specs{k,2}, S.zobs, 'o', 'MarkerSize', 5, ...
             'MarkerEdgeColor', cols(1,:), 'MarkerFaceColor', 'w', 'LineWidth', 1);
    end
    xline(0, 'k:');
    set(gca, 'YDir','reverse'); ylim(zl); box on; grid on
    xlabel(specs{k,3});
    if k == 1, ylabel('depth (m)'); else, set(gca,'YTickLabel',[]); end
    title(specs{k,4}, 'FontWeight','normal')
end

% ---- 3: p', the reconstructed profile --------------------------------
ax(3) = subplot(1,5,3); hold on
plot(S.pp, S.z, 'k-', 'LineWidth', 1.8);
xline(0, 'k:');
set(gca, 'YDir','reverse'); ylim(zl); box on; grid on
set(gca,'YTickLabel',[])
xlabel('p'' (Pa)')
title('perturbation pressure', 'FontWeight','normal')

% ---- 4: depth-resolved flux density ----------------------------------
% Faint: this hour. Heavy: the phase average, which is the one that integrates
% to the black arrow in panel 5.
ax(4) = subplot(1,5,4); hold on
fade = @(c) 1 - 0.6*(1 - c);        % toward white; plainer than an alpha channel
hi1 = plot(S.upp_inst, S.z, '-', 'LineWidth',0.9, 'Color',fade(cols(1,:)));
hi2 = plot(S.vpp_inst, S.z, '-', 'LineWidth',0.9, 'Color',fade(cols(2,:)));
ha1 = plot(S.upp,      S.z, '-', 'LineWidth',2.2, 'Color',cols(1,:));
ha2 = plot(S.vpp,      S.z, '-', 'LineWidth',2.2, 'Color',cols(2,:));
xline(0, 'k:');
set(gca, 'YDir','reverse'); ylim(zl); box on; grid on
set(gca,'YTickLabel',[])
xlabel('flux density (W m^{-2})')
legend([ha1 ha2 hi1 hi2], {'\langle u''p'' \rangle','\langle v''p'' \rangle', ...
       'u''p'', this hour','v''p'', this hour'}, 'Location','west', 'Box','off')
title(sprintf('flux density (%d-h mean)', numel(S.iwin)), 'FontWeight','normal')

% ---- 5: the depth-integrated flux -----------------------------------
ax(5) = subplot(1,5,5); hold on
all_mag = [hypot(S.F(1),S.F(2)), hypot(S.F_inst(1),S.F_inst(2)), ...
           hypot(S.F_mean(1),S.F_mean(2)), hypot(S.F_mode(:,1),S.F_mode(:,2))'];
lim = 1.2*max([all_mag, eps]);

qm = gobjects(1,nshow);
for j = 1:nshow
    qm(j) = quiver(0, 0, S.F_mode(j,1), S.F_mode(j,2), 0, 'LineWidth',1.0, ...
                   'MaxHeadSize',0.6, 'Color', cols(j+1,:));
end
qr = quiver(0, 0, S.F_mean(1), S.F_mean(2), 0, 'LineWidth',1.2, ...
            'MaxHeadSize',0.6, 'Color',[.55 .55 .55]);
qi = quiver(0, 0, S.F_inst(1), S.F_inst(2), 0, '--', 'LineWidth',1.2, ...
            'MaxHeadSize',0.6, 'Color',[.25 .25 .25]);
qp = quiver(0, 0, S.F(1), S.F(2), 0, 'k', 'LineWidth',2.4, 'MaxHeadSize',0.7);
xline(0,'k:'); yline(0,'k:');
axis equal; axis([-lim lim -lim lim]); box on; grid on
xlabel('F_x (W m^{-1})'); ylabel('F_y (W m^{-1})')
legend([qp qi qr qm], [{'F, phase-averaged','this hour, raw','record mean'} ...
        arrayfun(@(j)sprintf('mode %d',j), 1:nshow, 'uni',0)], ...
       'Location','southoutside', 'Box','off')
title(sprintf('|F| = %.2f kW/m toward %.0f\\circ', ...
      hypot(S.F(1),S.F(2))/1e3, out.flux.dir(it)), 'FontWeight','normal')

linkaxes(ax(1:4), 'y');
sgtitle(sprintf('%s   %s   %s band', out.name, datestr(S.time, 'dd mmm yyyy HH:MM'), ...
        cfg.band_name), 'FontWeight','bold')

S.ax = ax;
end

%% ======================================================================
function Xb = band_rows(X, cfg)
% Bandpass each row to cfg.M2_cpd, tolerating rows that are partly or wholly
% empty. bandpass_alford refuses NaNs outright, and an all-NaN sensor row
% survives fillmissing as all-NaN, so the two cases are separated here rather
% than left to error out of a plotting call.
Xb   = nan(size(X));
good = any(isfinite(X), 2);
if ~any(good), return; end
Xf = fillmissing(X(good,:), 'linear', 2, 'EndValues','nearest');
Xb(good,:) = bandpass_alford(Xf, cfg.dt_hr, cfg.M2_cpd, cfg.bw_c);
Xb(~isfinite(X)) = NaN;
end
