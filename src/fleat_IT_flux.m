function out = fleat_IT_flux(cfg, opt)
% FLEAT_IT_FLUX  One FLEAT mooring, from the gridded product to <u'p'> = F.
%
%   out = fleat_IT_flux(cfg)
%   out = fleat_IT_flux(cfg, opt)
%
%   cfg   from f1_config, f2_config, f3_config, f4_config or f6_config. The
%         config is the ONLY thing that distinguishes one mooring from another:
%         there is no mooring name anywhere in this function.
%
%   out   .name .lat .lon .H           where the arrow goes on the map
%         .time   [1 x nt]             hourly, datenum
%         .flux                        compute_tidal_flux output: Fx, Fy, mag,
%                                      dir, Fx_mode, Fy_mode, C, nflux [W/m]
%         .diag                        the numbers worth comparing across
%                                      moorings and the checks that should pass
%         .cfg                         the config as RUN, including any
%                                      overrides opt applied
%         .ok, .err                    false and a message if it failed
%         opt.keep = 'full' also attaches .Mh .bkgd .eta .modes .ahat .uv
%
%   out is scalar and every field is fixed in name and shape, so a loop over the
%   array can do R(i) = fleat_IT_flux(cfgs{i}) and get a struct array back.
%
%   ---------------------------------------------------------------------
%   OPTIONS
%   ---------------------------------------------------------------------
%   .verbose      print the running commentary. Default true.
%   .keep         'flux' (default) or 'full'. See the note on memory below.
%   .showPlots    the standard figures. Default false.
%   .figbase      add this to every figure number. Default 0. A loop over the
%                 array MUST set this (100*i) or mooring 5 draws over mooring 1.
%   .uv_nmodes    baroclinic velocity modes. Default cfg.nbc. This is the hook
%                 the degrees-of-freedom note in fit_uv_modes asks for: three
%                 ADCPs do not support six velocity modes, and the flux sum
%                 collapses to min(np, nu) so narrowing it here is safe.
%   .flux_periods phase-average window in tidal periods. Default cfg.flux_periods.
%   .flux_lowpass low-pass instead of a boxcar. Default false.
%
%   DIAGNOSTICS, all off by default. Each one answers a question about THIS
%   mooring and each costs real time, so a five-mooring run leaves them off and
%   you turn them on for one mooring at a time.
%   .glance       the T/S/u/v section overview (psectime_glance_mooring)
%   .check_order  also fit filter_first (Alford's order) and report the
%                 difference by mode. The pipeline itself always fits first.
%   .nmode_sweep  0, or an nmax: true m-mode refits of eta for m = 1..nmax
%   .eta_spectra  the eta amplitude spectra with continuum shoulders
%   .mode_movie   write the modes movie
%   .ellipse      the tidal ellipse figure
%   .synth_date   a datenum, or [] -- compare_synth_real_uvp at that date
%
%   ---------------------------------------------------------------------
%   WHY keep = 'flux' IS THE DEFAULT
%   ---------------------------------------------------------------------
%   The fullfield products are 34 to 184 MB each and the hourly struct plus the
%   modes is tens of MB per mooring. Five of those alive at once in a driver
%   loop is avoidable: the flux, the position and the diagnostics are all the
%   array-scale work needs, and they are kilobytes. Ask for 'full' when you are
%   looking at one mooring.
%
%   ---------------------------------------------------------------------
%   WHAT THIS IS, PHYSICALLY
%   ---------------------------------------------------------------------
%   p' from vertical mode fits of eta, the isopycnal displacement (Nash 2005's
%   xi: displacement of an isopycnal from its time-mean position, here against a
%   25-day running mean so mesoscale heave is removed rather than fitted). The
%   same fit gives the horizontal velocity modes. The flux is then the
%   phase-averaged product, summed over modes.
%
%   Adapted from turner_fleat2_IT_flux, which stays on disk as the F2 scratch
%   workspace. Methods from Alford 2003 and Nash 2005.
%
%   author TEJ, with help from Claude
%   15 October 2026

%% ==== 0. options and guards ==========================================
if nargin < 2 || isempty(opt), opt = struct(); end
opt = default(opt, 'verbose',      true);
opt = default(opt, 'keep',         'flux');
opt = default(opt, 'showPlots',    false);
opt = default(opt, 'figbase',      0);
opt = default(opt, 'glance',       false);
opt = default(opt, 'check_order',  false);
opt = default(opt, 'nmode_sweep',  0);
opt = default(opt, 'eta_spectra',  false);
opt = default(opt, 'mode_movie',   false);
opt = default(opt, 'ellipse',      false);
opt = default(opt, 'synth_date',   []);
opt = default(opt, 'uv_nmodes',    []);
opt = default(opt, 'flux_periods', []);
opt = default(opt, 'flux_lowpass', false);

if ~isstruct(cfg) || ~isfield(cfg,'mooring')
    error('fleat_IT_flux:nocfg', ...
        ['First argument must be a mooring config struct, e.g. ' ...
         'fleat_IT_flux(f2_config()).']);
end
tag = cfg.mooring;
V   = opt.verbose;

% THE BAND. Everything downstream -- bandpass_alford inside fit_uv_modes, the
% phase-average window in compute_tidal_flux -- reads cfg.M2_cpd. Copying
% band_cpd onto it here is what makes a K1 run a one-line change in the config
% instead of an edit in three functions. See the note in fleat_config_base.
if isfield(cfg,'band_cpd') && ~isempty(cfg.band_cpd)
    cfg.M2_cpd = cfg.band_cpd;
end
if ~isfield(cfg,'band_name'), cfg.band_name = 'M2'; end

if isempty(opt.flux_periods)
    if isfield(cfg,'flux_periods'), opt.flux_periods = cfg.flux_periods;
    else,                           opt.flux_periods = 2; end
end
if isempty(opt.uv_nmodes), opt.uv_nmodes = cfg.nbc; end

newfig = @(k) figure(opt.figbase + k);

say(V, tag, '================================================================\n');
say(V, tag, '%s  lat %.4f  lon %.4f  H %.0f m  band %s (%.4f cpd)  nbc %d\n', ...
    tag, cfg.lat, cfg.lon, cfg.H, cfg.band_name, cfg.M2_cpd, cfg.nbc);

% ---- fields that cannot be guessed ----------------------------------
% These two reach build_modes and f_background_v2 respectively, where a NaN
% produces either an error from somewhere deep or, worse, a profile. Refuse
% here, where the message can name the file to edit.
if ~isfinite(cfg.H)
    error('fleat_IT_flux:noH', ...
        ['cfg.H is NaN for %s. Run fleat_moorings() to read the recorded ' ...
         'bottom depth off the deployment file, then paste it into HREC in ' ...
         'fleat_config_base.'], tag);
end
if ~isfinite(cfg.mld)
    error('fleat_IT_flux:nomld', ...
        ['cfg.mld is NaN for %s. It is the shallowest sensor that carries ' ...
         'SALINITY, not the shallowest sensor -- set it in %s_config and see ' ...
         'the mixed-layer note in fleat_config_base for what it costs.'], ...
         tag, lower(tag));
end

% ---- the banner that stops a placeholder passing for a result -------
if isfield(cfg,'unverified') && ~isempty(cfg.unverified)
    warning('fleat_IT_flux:unverified', ...
        ['%s is running with UNVERIFIED config fields: %s.\n' ...
         'The pipeline does not care and will return a flux anyway. ' ...
         'See the to-establish list in %s_config.'], ...
         tag, strjoin(cfg.unverified, ', '), lower(tag));
end

out = struct('name', tag, 'lat', cfg.lat, 'lon', cfg.lon, 'H', cfg.H, ...
             'time', [], 'flux', [], 'diag', struct(), 'cfg', cfg, ...
             'ok', false, 'err', '');
D = struct();        % accumulates out.diag

%% ==== 1. load, trim, QC, hourly, density =============================
if ~isfile(cfg.datafile)
    error('fleat_IT_flux:nodata', 'No gridded product for %s at\n  %s', tag, cfg.datafile);
end
say(V, tag, 'loading %s\n', cfg.datafile);
L = load(cfg.datafile);
if ~isfield(L, cfg.datavar)
    error('fleat_IT_flux:novar', ...
        '%s holds %s, not the "%s" expected by cfg.datavar.', ...
        cfg.datafile, strjoin(fieldnames(L)', ', '), cfg.datavar);
end
F = L.(cfg.datavar);  clear L

% Trim BEFORE the hourly average, so the background mean, N2, the modes, eta
% and every fit all see the same window. No-op when both ends are [].
F = trim_mooring_record(F, cfg);

% Range-gate the RAW samples, before averaging: one bad 10-minute value inside
% an hourly mean contaminates the whole hour. No-op when the gates are [].
[F.T, F.S, D.qc] = qc_ts(F.T, F.S, F.depth, cfg);

Mh = mooring_hourly_mean(F, cfg);
clear F                                  % the big one; nothing below needs it
Mh = compute_density_mooring(Mh, cfg);

% Mixed layer depth, two ways. This does NOT feed the modes -- f_background_v2
% still caps N^2 above cfg.mld -- it is the only observational check on that
% cap, and at these moorings it mostly reports that the layer is above the
% shallowest density sensor and so cannot be checked. Read the header.
Mh = compute_MLD_mooring(Mh, cfg);

out.time   = Mh.time(:).';
D.nt       = numel(out.time);
D.record   = [out.time(1) out.time(end)];
D.ndays    = out.time(end) - out.time(1);
D.nsens    = size(Mh.depth, 1);
say(V, tag, '%d hourly steps, %.0f days, %d sensors\n', D.nt, D.ndays, D.nsens);

% The MLD numbers worth carrying across moorings. The resolved FRACTION is the
% one to read first: a small one means the medians beside it describe the few
% hours the mixed layer was deep enough to see, not the record.
D.mld_zref       = Mh.ml.info.zref_med;
D.mld_res_frac   = Mh.ml.info.frac_res_t;
D.mld_thresh     = Mh.ml.info.z_thresh_med;
D.mld_pea        = Mh.ml.info.z_pea_med;
D.mld_pea_diff   = Mh.ml.info.pea_minus_thresh;
% The energy the threshold MLD implies (Reichl et al. 2022) and the PEA over a
% fixed 100 m below the reference. The second one is defined in every hour with
% a profile, so it is the only one of these safe to plot against time.
D.mld_phi_thresh = Mh.ml.info.phi_thresh_med;
D.mld_phi_span   = Mh.ml.info.phi_span_med;
say(V, tag, 'MLD resolved in %.0f%% of hours; medians %.0f m (threshold) %.0f m (PEA)\n', ...
    100*D.mld_res_frac, D.mld_thresh, D.mld_pea);

if opt.showPlots
    % Both estimates against cfg.mld, with the unresolved hours marked. If the
    % grey dominates, the cap is not being checked by this record.
    res  = isfinite(Mh.ml.z_thresh) & ~Mh.ml.unresolved;
    resp = isfinite(Mh.ml.z_pea) & ~Mh.ml.unresolved & ~Mh.ml.pea_pinned;
    newfig(7); set(gcf,'Color','w')
    subplot(1,3,[1 2])
      % Explicit handles into legend: the yline below is a ConstantLine and
      % would otherwise take a label slot and shift all four.
      h1 = plot(Mh.time, Mh.ml.zref, 'Color',[.6 .6 .6]); hold on; box on
      h2 = plot(Mh.time(~res), Mh.ml.z_thresh(~res), '.', 'Color',[.8 .8 .8], 'MarkerSize',4);
      h3 = plot(Mh.time(res),  Mh.ml.z_thresh(res),  '.', 'MarkerSize',6);
      h4 = plot(Mh.time(resp), Mh.ml.z_pea(resp),    '.', 'MarkerSize',6);
      if isfinite(cfg.mld), yline(cfg.mld, 'r--', 'cfg.mld'); end
      datetick('x','mmm'); set(gca,'YDir','reverse'); ylabel('z [m]')
      % plot() of an all-false selection returns an EMPTY handle, not an empty
      % line, so [h1 h2 h3 h4] silently loses entries and legend warns
      % "Ignoring extra legend entries" with the labels then attached to the
      % wrong curves. At F6 nothing resolves, so h3 and h4 are both empty.
      % Pair each handle with its label and keep only the pairs that exist.
      hh  = [h1 h2 h3 h4];
      lab = {'shallowest \sigma','unresolved (bound)','threshold','PEA'};
      kp  = ~cellfun(@isempty, {h1, h2, h3, h4});
      legend(hh, lab(kp), ...
             'Location','southoutside','Orientation','horizontal','Box','off')
      title(sprintf('%s mixed layer depth, resolved %.0f%% of hours', ...
                    tag, 100*D.mld_res_frac), 'FontWeight','normal')
    subplot(1,3,3)
      histogram(Mh.ml.z_thresh(res),  'Orientation','horizontal'); hold on; box on
      histogram(Mh.ml.z_pea(resp),    'Orientation','horizontal')
      set(gca,'YDir','reverse'); xlabel('hours'); title('resolved only','FontWeight','normal')
    linkaxes(findall(gcf,'Type','axes'),'y')
end

if opt.glance && opt.showPlots
    % The 7th argument is zlim (a max depth), NOT a showPlots flag -- the
    % signature changed. opt.figbase goes in 8th so an array run does not draw
    % every mooring into figures 1 and 2.
    psectime_glance_mooring(Mh, cfg, [], 34.4:.02:34.9, -.25:.01:.25, [], [], opt.figbase);
end

%% ==== 2. background density and stratification =======================
% f_background_v2 builds BOTH N2 estimates every run; cfg.N2_source selects
% which one lands in bkgd.N2full and so reaches build_modes. Flip it afterwards
% with set_N2_source -- both profiles are already on bkgd, so the switch costs
% nothing and does not rebuild the background. Everything from cell 3 on is
% stale after a flip.
Mh.bkgd = f_background_v2(Mh, cfg);
if opt.showPlots
    plot_background(Mh.bkgd, cfg, opt.figbase + 4);
    set(gcf, 'Position', [80 80 1100 420])
end

% bkgd.N2full
%% ==== 3. isopycnal displacement ======================================
eta = f2_eta(Mh.depth, Mh.sigma, Mh.bkgd, Mh.time, cfg);

% Drop reference levels resolved too rarely to mean anything. bkgd and eta are
% cut together so the design matrix and the data keep the same rows.
keep = eta.cover(:) >= cfg.min_coverage;
if ~all(keep)
    say(V, tag, 'dropping %d level(s) with coverage < %.2f\n', nnz(~keep), cfg.min_coverage);
    Mh.bkgd.zref = Mh.bkgd.zref(keep);   Mh.bkgd.sig_t = Mh.bkgd.sig_t(keep);
    eta.eta  = eta.eta(keep,:);   eta.ziso  = eta.ziso(keep,:);
    eta.zbar = eta.zbar(keep,:);  eta.cover = eta.cover(keep);
end
D.nlevels = numel(Mh.bkgd.zref);

Mh.iso = struct('time', Mh.time, 'depth', eta.ziso, 'eta', eta.eta, 'zbar', eta.zbar);
Mh.iso.sig_t = Mh.bkgd.sig_t;
Mh.iso.cover = eta.cover;

% Time-mean eta should be ~0 by construction (the running mean defines the
% reference depth). It is a check, not a result: 0.31 m at F2.
D.eta_mean_max = max(abs(mean(eta.eta, 2, 'omitnan')));
D.eta_rms      = std(eta.eta, 0, 2, 'omitnan')';
say(V, tag, 'max |time-mean eta| %.2f m (want ~0); rms eta %.1f to %.1f m\n', ...
    D.eta_mean_max, min(D.eta_rms), max(D.eta_rms));

if opt.showPlots
    % 10 days: individual tides visible
    i0 = max(1, round(D.nt/3));
    newfig(3); set(gcf,'Color','w')
    psecline_mooring(Mh.iso, 'eta', [], 'depth', [], cfg, i0:min(D.nt, i0+10*24));
end

%% ==== 4. modes on the reference grid =================================
M = build_modes(Mh.bkgd.zfull, Mh.bkgd.N2full, cfg.lat, cfg.nbc, cfg.rho0);
D.c = M.c(:)';
say(V, tag, 'eigenspeeds: %s m/s\n', mat2str(round(M.c(:)'*1000)/1000));

if opt.showPlots
    newfig(5)
    subplot(1,3,1)
      plot(M.Zeta, M.z, 'LineWidth',1.4); hold on
      plot(zeros(size(Mh.bkgd.zref)), Mh.bkgd.zref, 'k.', 'MarkerSize',9)
      xlabel('\zeta modes'); ylabel('z (m)'); title('displacement (for \eta)')
      legend(arrayfun(@(j)sprintf('%d',j),1:cfg.nbc,'uni',0),'Location','southeast','Box','off')
      set(gca,'YDir','reverse')
    subplot(1,3,2)
      plot(M.Zu, M.z, 'LineWidth',1.4)
      xlabel('u modes'); title('horizontal velocity (for u'')')
      set(gca,'YDir','reverse')
    subplot(1,3,3)
      plot(M.P, M.z, 'LineWidth',1.4)
      xlabel('P'); title('pressure')
      set(gca,'YDir','reverse')
    sgtitle(sprintf('%s modes', tag))
end

%% ==== 5. fit and bandpass ============================================
fopt = struct('nmodes',cfg.nbc, 'weight',cfg.fit_weight, 'minobs',cfg.fit_minobs, ...
              'maxcond',cfg.fit_maxcond);

% FIT FIRST, ALWAYS. The modal fit is linear in eta at each timestep, so
% bandpassing eta and then fitting is algebraically the same operation as
% fitting raw eta and then bandpassing the amplitude series. The two orders
% commute, which means the choice between them is not physics -- it is only a
% question of where the GAPS have to be filled, because filtfilt is recursive
% and cannot take a NaN (bandpass_alford.m:36).
%
% Fitting first puts the fill on cfg.nbc smooth amplitude series. Filtering
% first puts it on every sensor level, which at F2 is 40 gappy rows. So fitting
% first is right at every mooring in this array, cfg.band_order was removed from
% fleat_config_base on 10 Nov 2026, and the guard below refuses a stale one
% rather than letting it run silently.
%
% Alford's own order (filter at the sensor levels, then fit) still exists as the
% opt.check_order diagnostic below, which runs it and reports the relative
% difference by mode -- 8-14% at F2. That is the way to trigger it. If you ever
% need filter_first as the PRODUCTION path again, it is four lines: eta_M2 below
% is already the band-passed eta, so refit from it. The full version is kept
% commented in turner_fleat2_IT_flux.m, cell 5.
if isfield(cfg,'band_order') && ~isempty(cfg.band_order) && ...
   ~strcmp(cfg.band_order, 'fit_first')
    error('fleat_IT_flux:band_order', ...
        ['cfg.band_order = "%s" is no longer honoured: this pipeline always ' ...
         'fits first (see the note above, and opt.check_order for the ' ...
         'comparison). Remove the field from the config.'], cfg.band_order);
end

% Fit raw eta (NaN-tolerant), then filter the amplitude series. Gap-filling
% happens on nbc series instead of on every sensor level.
[ahat_raw, fit] = fit_eta_modes(eta.eta, Mh.bkgd.zref, M, fopt);
gap  = ~all(isfinite(ahat_raw),1);
afil = fillmissing(ahat_raw, 'linear', 2, 'EndValues','nearest');
ahat_M2 = bandpass_alford(afil, cfg.dt_hr, cfg.M2_cpd, cfg.bw_c);
ahat_f  = bandpass_alford(afil, cfg.dt_hr, cfg.f_cpd,  cfg.bw_c);
% Re-mask the hours the fill invented, so no interpolated hour is ever reported
% as a result. The filter has already run; nothing downstream re-filters these.
ahat_M2(:,gap) = NaN;  ahat_f(:,gap) = NaN;

% Banded eta at the reference levels, for the sweep in cell 6 and for anything
% comparing the fit against the observations. Not an input to the fit.
eta_fill = fillmissing(eta.eta, 'linear', 2, 'EndValues','nearest');
eta_M2   = bandpass_alford(eta_fill, cfg.dt_hr, cfg.M2_cpd, cfg.bw_c);
eta_M2(~isfinite(eta.eta)) = NaN;

D.nfit_complete = nnz(all(isfinite(ahat_M2),1));
D.rms_band      = std(ahat_M2, 0, 2, 'omitnan')';
D.rms_f         = std(ahat_f,  0, 2, 'omitnan')';
say(V, tag, 'complete fits at %d of %d timesteps (%.0f%%)\n', ...
    D.nfit_complete, D.nt, 100*D.nfit_complete/D.nt);
say(V, tag, 'rms %s amplitude by mode:   %s m\n', cfg.band_name, mat2str(round(D.rms_band,2)));
say(V, tag, 'rms near-f amplitude by mode: %s m\n', mat2str(round(D.rms_f,2)));

% ---- does the filter/fit order matter here? -------------------------
% THE filter_first PATH. Alford's order: filter at the sensor levels, then fit
% each timestep. Reported as a relative difference by mode against the fit_first
% result above. If the two agree to a percent or so, the choice is immaterial for
% this mooring and you can say so in one sentence instead of defending it. It is
% 8-14% at F2, where the fill lands on 40 gappy levels instead of 6 amplitudes.
%
% eta_M2 is already the banded eta, so the whole alternate order is one fit.
D.band_order_rel = [];
if opt.check_order
    ahat_alt = fit_eta_modes(eta_M2, Mh.bkgd.zref, M, fopt);
    D.band_order_rel = (std(ahat_M2 - ahat_alt, 0, 2, 'omitnan') ./ ...
                        std(ahat_M2, 0, 2, 'omitnan'))';
    say(V, tag, 'fit_first vs filter_first, rms relative difference by mode: %s\n', ...
        mat2str(round(D.band_order_rel*1000)/1000));
end

if opt.showPlots
    newfig(6)
    for j = 1:cfg.nbc
        subplot(cfg.nbc,1,j)
        plot(Mh.time, ahat_M2(j,:)); hold on
        c = orderedcolors("gem"); plot(Mh.time, ahat_f(j,:), 'Color', c(2,:))
        datetick('x','mmm'); ylabel(sprintf('eta m%d (m)',j)); box on
        if j==1
            legend(cfg.band_name,'near-f','Box','off','Orientation','horizontal')
            title(sprintf('%s modal displacement amplitudes', tag),'FontWeight','normal')
        end
    end
    linkaxes(findall(gcf,'Type','axes'),'y')
end

if opt.eta_spectra && opt.showPlots
    D.spectra = f2_eta_spectra(ahat_raw, Mh.time, cfg, ...
                               struct('shoulders',[1.25 1.60 2.55 3.50]));
end

%% ==== 6. is nbc the right number here? ===============================
% Two views of one question, per mooring. This is the field that most deserves
% to be checked separately at each site: nbc = 6 was chosen in 1515 m of water
% with 43 sensors, and F6 has 434 m and far fewer.
%
%   varexp     cumulative weighted variance of eta explained by modes 1..m.
%              Where it flattens, extra modes buy nothing.
%   residual   rms residual at each depth after subtracting modes 1..m. Still
%              structured at m = nbc means there is signal being missed.
D.varexp = median(fit.varexp, 2, 'omitnan')';
say(V, tag, 'median cumulative variance explained (vs the fitted field): %s\n', ...
    mat2str(round(D.varexp,3)));

if opt.showPlots
    newfig(8)
    subplot(1,2,1)
      plot(1:cfg.nbc, D.varexp, 'o-', 'LineWidth',1.4); grid on
      xlabel('modes retained'); ylabel('median cumulative variance explained')
      ylim([0 1]); xticks(1:cfg.nbc); title('where does this flatten?')
    subplot(1,2,2)
      hold on; box on
      for m = 1:cfg.nbc
          r = eta.eta - fit.A(:,1:m)*ahat_raw(1:m,:);
          plot(std(r, 0, 2, 'omitnan'), Mh.bkgd.zref, 'o-', 'LineWidth',1.2)
      end
      plot(std(eta.eta,0,2,'omitnan'), Mh.bkgd.zref, 'k--', 'LineWidth',1.2)
      xlabel('rms residual \eta [m]'); ylabel('z [m]')
      set(gca,'YDir','reverse')
      legend([arrayfun(@(m)sprintf('1-%d',m),1:cfg.nbc,'uni',0) {'no fit'}], ...
             'Location','northeast')
    sgtitle(sprintf('%s: how many modes?', tag))
end

% ---- true m-mode refits, raw and band ------------------------------
% The varexp above is cumulative within ONE nbc-mode fit. This refits for each
% m, which is the honest version of the question, and does it against raw eta
% and band-passed eta separately. Raw includes the broadband continuum that m
% modes were never going to represent, so the two numbers are very different
% and only the band one speaks to the flux. 0.665 vs ~0.86 at F2.
% ---------------------------------------------------------------------
% THE GATE BELONGS ON THIS CURVE. Added 10 Nov 2026, after F1.
%
% The sweep as originally written refits with a bare backslash and NO
% conditioning test, so its curve rises monotonically in m whatever the geometry
% can actually identify -- the same way R^2 improves when you add predictors to
% a regression. At F1 that read 0.884 at m=3 and 0.912 at m=4 while the real
% gated fit delivered 0.767 and printed NEGATIVE cumulative variance explained
% at m=1 and m=2, the signature of four large amplitudes that only cancel once
% they are all summed. Variance explained alone cannot see that.
%
% So the sweep now reports three curves against m, and the gap between the first
% two IS the conditioning:
%   pct      cumulative variance explained, NO gate  (what the sweep always did)
%   pct_cond the same refit, scoring only the hours that pass cfg.fit_maxcond
%   frac     the fraction of hours that pass -- the hidden cost of the second
%            curve, because 0.9 computed on 3% of hours is not a result
% plus cond(A_w) itself, median and 10-90%, which is the thing being traded.
%
% Ranges differ and that is the point: the first three live in (-Inf, 1] and the
% last in [1, Inf), so cond gets its own log axis in the figure. cond RISES with
% m -- strictly, since dropping a column cannot raise the largest singular value
% nor lower the smallest, so cond(A(:,1:m)) <= cond(A(:,1:m+1)) always. Read the
% crossing: the largest m whose cond still clears maxcond on most hours.
D.pct_nmodes = [];  D.pct_nmodes_cond = [];
D.cond_nmodes = []; D.frac_nmodes = [];
if opt.nmode_sweep > 0
    nmax = opt.nmode_sweep;
    Mx   = build_modes(Mh.bkgd.zfull, Mh.bkgd.N2full, cfg.lat, nmax, cfg.rho0);
    Z    = interp1(Mx.z, Mx.Zeta, Mh.bkgd.zref);     % modes at the fit levels
    w    = fit.w(:);                                 % same weights as the fit
    sets = {eta.eta, eta_M2};
    maxc = cfg.fit_maxcond;
    pct  = nan(nmax, 2);   pctc = nan(nmax, 2);
    cstat= nan(nmax, 3);   frac = nan(nmax, 1);
    % The gate is nmax+2 for every m, not m+2, so all m are scored on the SAME
    % hours and the curve is a comparison rather than a mix of subsamples. The
    % cost is that nmax alone decides whether the sweep runs at all: at F6,
    % nmode_sweep = 10 against 7 reference levels asked for 12 good levels, no
    % hour qualified, den stayed 0, and 1 - 0/0 printed ten NaNs for a four-mode
    % run. Report that instead of returning it.
    nhr = zeros(1,2);

    % ---- pass 1: conditioning per qualifying hour, per m -------------
    % cond depends only on WHICH ROWS survive and on m, never on the data
    % values, so one pass serves both sets -- see the equality check below.
    okmat = isfinite(eta.eta);
    if ~isequal(isfinite(eta_M2), okmat)
        error('fleat_IT_flux:sweep_mask', ...
            ['the sweep assumes eta_M2 is masked exactly where eta.eta is NaN ' ...
             '(cell 5) so that cond is shared between the two sets. That is no ' ...
             'longer true -- compute cond per set before trusting this curve.']);
    end
    iq   = find(sum(okmat,1) >= nmax + 2);
    cn_m = nan(numel(iq), nmax);
    for k = 1:numel(iq)
        ok = okmat(:,iq(k));  sw = sqrt(w(ok));
        for m = 1:nmax
            cn_m(k,m) = cond(sw .* Z(ok,1:m));
        end
    end
    if ~isempty(iq)
        cstat = [prctile(cn_m,10,1); median(cn_m,1); prctile(cn_m,90,1)]';
        frac  = mean(cn_m <= maxc, 1)';
    end

    % ---- pass 2: variance explained, with and without the gate -------
    for s = 1:2
        X = sets{s};
        num  = zeros(nmax,1);  den  = 0;
        numg = zeros(nmax,1);  deng = zeros(nmax,1);
        for k = 1:numel(iq)
            it = iq(k);
            e  = X(:,it);  ok = okmat(:,it);
            sw = sqrt(w(ok));
            v0 = sum(w(ok).*e(ok).^2);
            den = den + v0;
            nhr(s) = nhr(s) + 1;
            for m = 1:nmax
                a = (sw.*Z(ok,1:m)) \ (sw.*e(ok));
                r = sum(w(ok).*(e(ok) - Z(ok,1:m)*a).^2);
                num(m) = num(m) + r;
                if cn_m(k,m) <= maxc                 % the gated accounting
                    numg(m) = numg(m) + r;  deng(m) = deng(m) + v0;
                end
            end
        end
        if den > 0
            pct(:,s) = 1 - num/den;
        end
        g = deng > 0;                                % NaN where no hour passed
        pctc(g,s) = 1 - numg(g)./deng(g);
    end

    if all(nhr == 0)
        warning('fleat_IT_flux:sweep_starved', ...
            ['%s: opt.nmode_sweep = %d needs %d good levels in an hour and ' ...
             'zref has only %d, so no hour qualifies and the sweep is all ' ...
             'NaN. Lower it to %d or less.'], ...
             tag, nmax, nmax+2, D.nlevels, max(D.nlevels-2, 1));
    elseif min(nhr) < 0.1*D.nt
        say(V, tag, 'nmode_sweep scored only %d raw / %d band hours of %d (needs %d levels)\n', ...
            nhr(1), nhr(2), D.nt, nmax+2);
    end
    D.pct_nmodes      = pct;         % [nmax x 2], columns raw then band
    D.pct_nmodes_cond = pctc;        % the same, gated at cfg.fit_maxcond
    D.cond_nmodes     = cstat;       % [nmax x 3], p10 / median / p90
    D.frac_nmodes     = frac;        % [nmax x 1], hours passing the gate
    say(V, tag, '%% explained, raw:  %s\n', mat2str(round(pct(:,1)',3)));
    say(V, tag, '%% explained, %4s: %s\n', cfg.band_name, mat2str(round(pct(:,2)',3)));
    say(V, tag, '%% explained, %4s, gated at cond<=%g: %s\n', ...
        cfg.band_name, maxc, mat2str(round(pctc(:,2)',3)));
    say(V, tag, 'median cond by mode: %s\n', mat2str(round(cstat(:,2)')));
    say(V, tag, 'hours passing the gate: %s\n', mat2str(round(frac'*100)/100));

    if opt.showPlots
        Seta = struct('m', 1:nmax, 'pct', pct(:,2), 'pctc', pctc(:,2), ...
                      'frac', frac, 'cond', cstat, 'pct_extra', pct(:,1), ...
                      'pct_name',   [cfg.band_name '-band \eta varexp, no gate'], ...
                      'pctc_name',  [cfg.band_name '-band \eta varexp, gated'], ...
                      'extra_name', 'raw \eta varexp, no gate');
        pmode_sweep(Seta, maxc, ...
            sprintf('%s displacement: how many modes do these levels support?', tag), ...
            opt.figbase + 9);
    end

    % Carried so cell 9 can build the velocity sweep on the same mode set
    % without a second build_modes call.
    Mx_sweep = Mx;
end



%% ==== 8. the modes movie (optional) =================================
if opt.mode_movie
    % Per-timestep varexp against the BAND-passed eta, which is what the movie
    % should report; fit.varexp is against whatever the fit saw.
    fitM = fit;  fitM.varexp = nan(cfg.nbc, D.nt);
    for it = 1:D.nt
        e = eta_M2(:,it); ok2 = isfinite(e);
        if nnz(ok2) < cfg.fit_minobs || ~all(isfinite(ahat_M2(:,it))), continue; end
        v0 = sum(fit.w(ok2).*e(ok2).^2);
        for m = 1:cfg.nbc
            r = e(ok2) - fit.A(ok2,1:m)*ahat_M2(1:m,it);
            fitM.varexp(m,it) = 1 - sum(fit.w(ok2).*r.^2)/v0;
        end
    end
    mopt = struct('outfile', sprintf('modes_%s_%s.mp4', tag, cfg.band_name), ...
                  'fps', 8, 'pplot', true);
    % 200 days in is clear of the filter edges at F2, but F4 only has 164 days
    % of record, so fall back to two weeks past the first complete fit.
    first = find(all(isfinite(ahat_M2),1), 1);
    i0 = first + 200*24;
    if isempty(i0) || i0 + 10*24 > D.nt, i0 = first + 14*24; end
    if isempty(i0) || i0 + 2 > D.nt
        warning('fleat_IT_flux:nomovie', ...
            '%s has no 10-day window of complete fits to animate.', tag);
        mopt.frames = [];
    else
        mopt.frames = i0:2:min(D.nt, i0 + 10*24);       % 10 days, every 2 h
    end
    if ~isempty(mopt.frames)
        mode_movie(eta_M2, Mh.bkgd.zref, M, ahat_M2, fitM, Mh.time, mopt);
    end
end

%% ==== 9. velocity modes ==============================================
% Same baroclinic set as the displacement fit, plus mode 0. opt.uv_nmodes is
% the place to narrow the velocity fit if the ADCP count says to -- see the
% degrees-of-freedom note in fit_uv_modes. Row 1 of uhat/vhat is BAROTROPIC.
% minobs is passed explicitly: fit_uv_modes defaults it to 8 and does NOT read
% cfg.fit_minobs, so a mooring with fewer than 8 velocity levels (F6 has 7)
% rejected every hour silently while the eta side, which does get cfg, ran.
%
% maxcond is passed for the SAME REASON, added 10 Nov 2026. fit_uv_modes
% defaults it to 100 and does not read cfg either, so cfg.fit_maxcond silently
% governed only half the pipeline: an F1 run with the eta gate opened to Inf
% still had every velocity hour rejected at 100, returned an all-NaN uhat, and
% died in bandpass_alford three steps later. Tightening the config had the
% mirror-image problem -- the velocity fit stayed at 100 regardless.
UV = fit_uv_modes(Mh, M, cfg, struct('nmodes', opt.uv_nmodes, ...
                                     'minobs', cfg.fit_minobs, ...
                                     'maxcond', cfg.fit_maxcond));
D.nuv = size(UV.uhat_M2,1) - 1;
D.cond_uv = [prctile(UV.cond,10) median(UV.cond,'omitnan') prctile(UV.cond,90)];
say(V, tag, 'velocity cond at %d mode(s): median %.0f, 10-90%% %.0f-%.0f\n', ...
    opt.uv_nmodes, D.cond_uv(2), D.cond_uv(1), D.cond_uv(3));

% ---- the same three-curve sweep, for the VELOCITY fit ---------------
% Added 10 Nov 2026. The displacement sweep in cell 6 answers only half the
% question: at F1 the eta side cleared cond ~25 at m = 4 while the ADCP side
% rejected every hour at the same gate, so the binding constraint was the
% velocity matrix and nothing plotted it. opt.uv_nmodes is the knob this
% figure sets; compute_tidal_flux already sums min(np, nu) modes and says
% which one has no partner, so nbc and uv_nmodes may legitimately differ.
D.uv_sweep = [];
if opt.nmode_sweep > 0 && exist('Mx_sweep','var')
    D.uv_sweep = uv_mode_sweep(Mh, Mx_sweep, ...
        struct('nmax', opt.nmode_sweep, 'w', UV.w, ...
               'maxcond', cfg.fit_maxcond, 'minobs', cfg.fit_minobs));
    say(V, tag, 'uv %% explained: %s\n', mat2str(round(D.uv_sweep.pct',3)));
    say(V, tag, 'uv median cond by mode: %s\n', mat2str(round(D.uv_sweep.cond(:,2)')));
    say(V, tag, 'uv hours passing the gate: %s\n', ...
        mat2str(round(D.uv_sweep.frac'*100)/100));
    if D.uv_sweep.nhr == 0
        warning('fleat_IT_flux:uv_sweep_starved', ...
            ['%s: the velocity sweep needs %d good bins in an hour and no hour ' ...
             'has them. Lower opt.nmode_sweep to read it.'], ...
             tag, max(cfg.fit_minobs, opt.nmode_sweep+3));
    elseif opt.showPlots
        S = D.uv_sweep;
        % UNFILTERED, and say so. The displacement figure's main curve is
        % band-passed eta; this one is the mean-removed but unfiltered u,v,
        % because that is the field fit_uv_modes actually fits (it bandpasses
        % the AMPLITUDES afterwards, fit_first). So this curve is the analogue
        % of the GREY raw-eta curve on the displacement figure, not of the blue
        % band curve, and the two blue lines are not comparable in level.
        S.pct_name  = 'u,v varexp, unfiltered, no gate';
        S.pctc_name = 'u,v varexp, unfiltered, gated';
        pmode_sweep(S, cfg.fit_maxcond, ...
            sprintf('%s velocity: how many modes do these bins support?', tag), ...
            opt.figbase + 10);
    end
end

if opt.showPlots
    nrow = size(UV.uhat_M2, 1);
    c = orderedcolors("gem");
    newfig(12)
    for r = 1:nrow
        subplot(nrow, 1, r)
        plot(Mh.time, UV.uhat_M2(r,:), 'Color', c(1,:)); hold on
        plot(Mh.time, UV.vhat_M2(r,:), 'Color', c(2,:))
        datetick('x','mmm'); box on
        if r == 1
            ylabel('mode 0 (BT) [m/s]')
            legend('u','v','Box','off','Orientation','horizontal')
            title(sprintf('%s velocity amplitudes', tag),'FontWeight','normal')
        else
            ylabel(sprintf('mode %d [m/s]', r-1))
        end
    end
    linkaxes(findall(gcf,'Type','axes'),'y')
end

%% ==== 10. the flux ===================================================
% Record the window that actually ran, not the one the config file suggests.
% out.cfg is documented as the config AS RUN, and opt.flux_periods can override
% cfg.flux_periods -- anything reproducing the phase average afterwards
% (puvp_slice's depth-resolved <u'p'>) has to average over the same hours or it
% will not integrate back to out.flux.
cfg.flux_periods = opt.flux_periods;
cfg.flux_lowpass = opt.flux_lowpass;

out.flux = compute_tidal_flux(ahat_M2, UV.uhat_M2, UV.vhat_M2, M, cfg, ...
               struct('periods', opt.flux_periods, 'lowpass', opt.flux_lowpass));

D.Fx_mean  = mean(out.flux.Fx,  'omitnan');
D.Fy_mean  = mean(out.flux.Fy,  'omitnan');
D.mag_mean = mean(out.flux.mag, 'omitnan');
D.dir_mean = mod(90 - atan2d(D.Fy_mean, D.Fx_mean), 360);
% A record-mean much smaller than the typical magnitude means the flux reverses
% direction, i.e. a standing rather than a propagating wave. This ratio is the
% single most comparable number across the array.
D.steadiness = hypot(D.Fx_mean, D.Fy_mean) / D.mag_mean;
say(V, tag, 'steadiness |mean F| / mean |F| = %.2f (1 = purely propagating)\n', D.steadiness);

%% ==== 11. optional per-mooring diagnostics ===========================
if opt.ellipse && opt.showPlots
    D.ellipse = plot_tidal_ellipse(Mh, cfg);
end
% if ~isempty(opt.synth_date) && opt.showPlots
%     % synthetic (NW arrival) vs real, same modes
%     compare_synth_real_uvp(Mh, ahat_M2, cfg, opt.synth_date, "to_deg", 360);
% end

%% ==== 12. assemble ===================================================
out.diag = D;
out.cfg  = cfg;
out.ok   = true;

if strcmp(opt.keep, 'full')
    out.Mh    = Mh;
    out.bkgd  = Mh.bkgd;
    out.eta   = eta;
    % eta in the band, at the reference levels. Cell 5 builds this and would
    % otherwise throw it away; keeping it means anything comparing the fit
    % against the observations (puvp_slice) does not have to re-derive it and
    % risk filtering it differently from the way the fit saw it.
    out.eta.band = eta_M2;
    out.modes = M;
    out.ahat  = struct('band', ahat_M2, 'near_f', ahat_f, 'raw', ahat_raw, 'fit', fit);
    out.uv    = UV;
    % p' on the full mode grid, in the band. M.P is already in Pa per metre of
    % displacement amplitude (build_modes), so this is the whole calculation --
    % puvp_slice.m:32 states the same identity and needs the field at lines
    % 199 and 269-270 to form <u'p'>(z).
    %
    % The name this line carried until 30 Oct 2026 was pprime_A, which it
    % inherited from turner_fleat2_IT_flux_verbose.m:323 and which was never
    % assigned here -- so keep = 'full' errored at every mooring. In that
    % script "A" distinguished the modal pressure from a route B that
    % integrated the density perturbation hydrostatically,
    % g*cumtrapz(z, rho_p) minus its depth mean, as a check that M.P and the
    % hydrostatic integral agree. That check is NOT reinstated here and the
    % suffix is dropped with it: it is a verification of build_modes, not a
    % per-run diagnostic, and it has nothing to do with cfg.band_order.
    out.pprime = M.P * ahat_M2;
end

say(V, tag, 'done.\n');
end

%% ======================================================================
function opt = default(opt, f, v)
% Fill one option if the caller left it out. A caller who passes [] on purpose
% keeps the [], which is how uv_nmodes and flux_periods signal "take it from
% the config" a few lines further down.
if ~isfield(opt, f), opt.(f) = v; end
end

%% ======================================================================
function say(on, tag, fmt, varargin)
% Every line of commentary carries the mooring name. Without this, the log of a
% five-mooring loop is unreadable: compute_tidal_flux and f_background_v2 both
% print freely and nothing says which site they are talking about.
if on
    fprintf('[%s] ', tag);
    fprintf(fmt, varargin{:});
end
end
