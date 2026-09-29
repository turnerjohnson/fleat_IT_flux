% turner_fleat2_IT_flux.m

% Main workspace for calculating <u'p'> = F for Palau FLEAT moorings.

% First getting p' from vertical mode fits of η, isopycnal displacement, to
% get full-column values for density. Then integrate eta/density to get p.





% contains functions:
%    f2_config.m        sets values and assumptions for F2 mooring, incl
%                           site info, sampling info, bandpass filtering 
%                           preferences, N2 stratification smoothing info, 
%                           and mode fitting info.
%    f2_background.m    calculates mean background density profile and N2, 
%                           sets mixed-layer extrapolation for a full depth
%                           grid
%    f2_eta.m           calculates isopycnal displacement from moored T/S
%                           array, or Nash2005's \xi: vertical displacement
%                           of isopycnal relative to time-mean position.
%                           this uses running mean instead of recird-mean
%                           to remove mesoscale heave (~25d)


% adapted from turner_tinker_fleat.m

% author TEJ, based on methods from Alford03 and Nash05 with help from Claude
% 5 Sept 2026

% gridded T, S, z data from Anna Simpson
% gridded u, v, data from Conrad Luecke

% --------------------------------------------
clear; close all;
addpath(genpath('/Users/tej/Documents'))

%% 0.0 config and load data.
cfg = f2_config();
F2 = load("/Users/tej/Documents/research_fleat/data/simpson/F2_fullfield.mat").fullfield;

%% 0.0 alternate: sthentic dataset to check physics logic
[F2, cfg] = synth_baroclinic();
%% 0. take hourly average

k = cfg.avg_k; % k=6 samples per hour.
hravg = @(X) reshape( squeeze(mean(reshape( X(:,1:floor(size(X,2)/k)*k), ...
              size(X,1), k, []), 2, 'omitnan')), size(X,1), []);

F2_h.time = hravg(F2.time);
F2_h.depth = hravg(F2.depth);          % [nsens x nt] instrument depth, positive down
F2_h.T = hravg(F2.T);
F2_h.S = hravg(F2.S);    
F2_h.u = hravg(F2.u);   
F2_h.v = hravg(F2.v);   

fprintf('%d sensors, %d hours (%.0f days)\n', size(F2_h.depth,1), numel(F2_h.time), F2_h.time(end)-F2_h.time(1));


%% 1. compute density

F2_h.p = sw_pres(F2_h.depth, cfg.lat);
F2_h.sigma = sw_pden(F2_h.S, F2_h.T, F2_h.p, 0) - 1000;

% How much of the density signal is salinity doing? At 8.5 N with a barrier
% layer this is not rhetorical, and it is worth having the number ready if
% anyone asks how far you trust the glider S.

sigma_fixedS = sw_pden(repmat(mean(F2_h.S,2,'omitnan'),1,size(F2_h.S,2)), F2_h.T, F2_h.p, 0) - 1000;
fprintf('rms sigma difference, glider S vs time-mean S: %.4f kg m^-3\n', ...
        std(F2_h.sigma(:) - sigma_fixedS(:), 'omitnan'));

%% optional: look at section data

temp_cont = [];%2:1:25;
sal_cont = 34.4:.02:34.9;
uv_cont = -.25:.01:.25;
time_lim = 2900:2999;

% figure(1)
% subplot(4,1,1)
% psectime_mooring(F2_h, 'T', [], temp_cont, [], [], 'depth', [], cfg);
% subplot(4,1,2)
% %psectime_mooring(F2_h, 'S', [], sal_cont, [], [], 'depth',[], cfg);
% subplot(4,1,3)
% psectime_mooring(F2_h, 'u', [], uv_cont, [], [], 'depth', [], cfg);
% subplot(4,1,4)
% psectime_mooring(F2_h, 'v', [], uv_cont, [], [], 'depth', [], cfg);

figure(2)
psecline_mooring(F2_h, 'sigma',  [],  'depth', [], cfg);


figure(1)
subplot(4,1,1)
psectime_mooring(F2_h, 'T', [], temp_cont, [], [], 'depth', 300, cfg, time_lim);
subplot(4,1,2)
psectime_mooring(F2_h, 'S', [], sal_cont, [], [], 'depth', 300, cfg, time_lim);
subplot(4,1,3)
psectime_mooring(F2_h, 'u', [], uv_cont, [], [], 'depth', 300, cfg, time_lim);
subplot(4,1,4)
psectime_mooring(F2_h, 'v', [], uv_cont, [], [], 'depth', 300, cfg, time_lim);

% figure(2)
% subplot(4,1,1)
% psecline_mooring(F2_h, 'T', temp_cont, 'depth', 150, cfg, time_lim);
% subplot(4,1,2)
% psecline_mooring(F2_h, 'S',  sal_cont, 'depth', 150, cfg, time_lim);
% subplot(4,1,3)
% psecline_mooring(F2_h, 'u',  uv_cont, 'depth', 150, cfg, time_lim);
% subplot(4,1,4)
% psecline_mooring(F2_h, 'v',  uv_cont, 'depth', 150, cfg, time_lim);

%% optional: look at glider data.
% gli = load_spray_fleat("/Users/tej/Documents/research_fleat/data/spray_FLEAT.nc");
% 
% % just what's near F2 during the deployment
% gli = load_spray_fleat("/Users/tej/Documents/research_fleat/data/spray_FLEAT.nc", ...
%         'lat0', 8.5027, 'lon0', 134.5918, 'radius', 100, ...
%         'tlim', [F2_h.time(1) F2_h.time(end)]);
% 
% plot_glider_coverage(gli, 8.5027, 134.5918, 100, ...
%                      [F2_h.time(1) F2_h.time(end)], 'spray_FLEAT.nc')

%% 2. background/mean density, stratification profile
F2_h.bkgd = f2_background(F2_h.depth, F2_h.sigma, cfg);

figure('Position',[80 80 1100 420],'Color','w')
subplot(1,3,1)
     % mean sigma profile 
  plot(F2_h.bkgd.sig_full, F2_h.bkgd.zfull, 'LineWidth',1.3); hold on
     % stabilized mean sigma scattered on surfaces
  plot(F2_h.bkgd.sig_t, F2_h.bkgd.zref, 'o', 'MarkerSize',4)
  yline(cfg.mld,'r--','MLD'); xlabel('\sigma_\theta'); ylabel('z [m]')
  title('mean profile + mixed layer'); legend('full','binned','Location','southwest')
  set(gca, 'YDir', 'reverse');

subplot(1,3,2)
     % smoothed, positive N2 from sig_full
  semilogx(F2_h.bkgd.N2full, F2_h.bkgd.zfull, 'LineWidth',1.3); hold on
     % median N2 from gridded Simpson fullfield
  semilogx(median(F2.N2,2,'omitnan'), mean(F2_h.depth,2,'omitnan'), '.')
  xlabel('N^2 [s^{-2}]'); title('N^2'); legend('from sigma','F.N2 median','Location','southeast')
  set(gca, 'YDir', 'reverse');

subplot(1,3,3)
     % buoyancy period: 60 * 2π / N2
  plot(2*pi./sqrt(F2_h.bkgd.N2full)/60, F2_h.bkgd.zfull, 'LineWidth',1.3)
  xlabel('buoyancy period [min]'); title('sanity check')
  set(gca, 'YDir', 'reverse');

% The two N2 curves should agree. 

%% 3. displacement of isopycnal surfaces
eta = f2_eta(F2_h.depth, F2_h.sigma, F2_h.bkgd, F2_h.time, cfg);

keep = eta.cover(:) >= cfg.min_coverage;
if ~all(keep)
    fprintf('dropping %d level(s) with coverage < %.2f\n', nnz(~keep), cfg.min_coverage);
    F2_h.bkgd.zref = F2_h.bkgd.zref(keep);  F2_h.bkgd.sig_t = F2_h.bkgd.sig_t(keep);
    eta.eta  = eta.eta(keep,:);   eta.ziso  = eta.ziso(keep,:);
    eta.zbar = eta.zbar(keep,:);  eta.cover = eta.cover(keep);
end

F2_h.iso = struct('time', F2_h.time, 'depth', eta.ziso, 'eta', eta.eta, 'zbar', eta.zbar);
F2_h.iso.sig_t = F2_h.bkgd.sig_t;
F2_h.iso.cover = eta.cover;

% 10-day window: individual tides visible
i0 = find(F2_h.time >= datenum(2016,10,1), 1);
figure('Color','w')
psecline_mooring(F2_h.iso, 'eta', [], 'depth', [], cfg, i0:i0+10*24);

%%  4. modes on reference grid (2m)
F2_h.modes = build_modes(F2_h.bkgd.zfull, F2_h.bkgd.N2full, cfg.lat, cfg.nbc, cfg.rho0);
fprintf('eigenspeeds: %s m/s\n', mat2str(round(F2_h.modes.c*1000)/1000));

figure('Position',[80 80 1000 420],'Color','w')
subplot(1,3,1)
  plot(F2_h.modes.Zeta, F2_h.modes.z, 'LineWidth',1.4); hold on
  plot(zeros(size(F2_h.bkgd.zref)), F2_h.bkgd.zref, 'k.', 'MarkerSize',9)
  xlabel('\zeta modes'); ylabel('z (m)'); title('displacement (for \eta)')
  legend(arrayfun(@(j)sprintf('%d',j),1:cfg.nbc,'uni',0),'Location','southeast','Box','off')
  set(gca, 'YDir', 'reverse');
subplot(1,3,2)
  plot(F2_h.modes.Zu, F2_h.modes.z, 'LineWidth',1.4)
  xlabel('u modes'); title('horizontal velocity (for u'')')
  set(gca, 'YDir', 'reverse');
subplot(1,3,3)
  plot(F2_h.modes.P, F2_h.modes.z, 'LineWidth',1.4)
  xlabel('P '); title('pressure')
  set(gca, 'YDir', 'reverse');



%% 5. fit and bandpass 
opt = struct('nmodes',cfg.nbc, 'weight',cfg.fit_weight, 'minobs',cfg.fit_minobs);

switch cfg.band_order
  case 'fit_first'
    % Fit raw eta (NaN-tolerant), then filter the amplitude series.
    [ahat_raw, fit] = fit_eta_modes(eta.eta, F2_h.bkgd.zref, F2_h.modes, opt);
    gap  = ~all(isfinite(ahat_raw),1);
    afil = fillmissing(ahat_raw, 'linear', 2, 'EndValues','nearest');
    ahat_M2 = bandpass_alford(afil, cfg.dt_hr, cfg.M2_cpd, cfg.bw_c);
    ahat_f  = bandpass_alford(afil, cfg.dt_hr, cfg.f_cpd,  cfg.bw_c);
    ahat_M2(:,gap) = NaN;  ahat_f(:,gap) = NaN;

    eta_fill = fillmissing(eta.eta, 'linear', 2, 'EndValues','nearest');
    eta_M2   = bandpass_alford(eta_fill, cfg.dt_hr, cfg.M2_cpd, cfg.bw_c);
    eta_M2(~isfinite(eta.eta)) = NaN;

  case 'filter_first'
    % Alford's order: filter at the sensor levels, then fit each timestep.
    eta_fill = fillmissing(eta.eta, 'linear', 2, 'EndValues','nearest');
    eta_M2   = bandpass_alford(eta_fill, cfg.dt_hr, cfg.M2_cpd, cfg.bw_c);
    eta_f    = bandpass_alford(eta_fill, cfg.dt_hr, cfg.f_cpd,  cfg.bw_c);
    eta_M2(~isfinite(eta.eta)) = NaN;  eta_f(~isfinite(eta.eta)) = NaN;
    [ahat_M2, fit] = fit_eta_modes(eta_M2, F2_h.bkgd.zref, F2_h.modes, opt);
    ahat_f         = fit_eta_modes(eta_f,  F2_h.bkgd.zref, F2_h.modes, opt);
    ahat_raw       = fit_eta_modes(eta.eta,  F2_h.bkgd.zref, F2_h.modes, opt);  % for the varexp plot

  otherwise
    error('unknown cfg.band_order "%s"', cfg.band_order);
end

fprintf('complete fits at %d of %d timesteps\n', nnz(all(isfinite(ahat_M2),1)), numel(F2_h.time));
fprintf('rms M2 amplitude by mode:     %s m\n', mat2str(round(std(ahat_M2,0,2,'omitnan')',2)));
fprintf('rms near-f amplitude by mode: %s m\n', mat2str(round(std(ahat_f, 0,2,'omitnan')',2)));

% --- settle the ordering question empirically, once ---
% Run the other order and see whether it changes anything. If the two agree to a
% percent or so, the choice is immaterial for this dataset and you can say so in
% one sentence instead of defending it.
cfg2 = cfg;
cfg2.band_order = 'filter_first';
if strcmp(cfg.band_order,'fit_first')
    eta_f2 = fillmissing(eta.eta,'linear',2,'EndValues','nearest');
    eta_f2 = bandpass_alford(eta_f2, cfg.dt_hr, cfg.M2_cpd, cfg.bw_c);
    eta_f2(~isfinite(eta.eta)) = NaN;
    ahat_alt = fit_eta_modes(eta_f2, F2_h.bkgd.zref, F2_h.modes, opt);
    rel = std(ahat_M2 - ahat_alt, 0, 2, 'omitnan') ./ std(ahat_M2, 0, 2, 'omitnan');
    fprintf('fit_first vs filter_first, rms relative difference by mode: %s\n', ...
            mat2str(round(rel'*1000)/1000));
end

figure('Position',[80 80 1000 550],'Color','w')
for j = 1:cfg.nbc
    subplot(cfg.nbc,1,j)
    %plot(F2_h.time, ahat_M2(j,:)); hold on; 
    c = orderedcolors("gem"); plot(F2_h.time, ahat_f(j,:), 'Color', c(2,:))
    datetick('x','mmm'); ylabel(sprintf('eta m%d (m)',j)); box on
    if j==1, legend('M_2','near-f','Box','off','Orientation','horizontal'); end
    linkaxes(findall(gcf,'Type','axes'),'y')
end
%% spectra sanity check
S = f2_eta_spectra(ahat_raw, F2_h.time, cfg, struct('shoulders',[1.25 1.60 2.55 3.50]));

%% ==== 6b. is nbc the right number? ==================================
% Data-driven, no simulation. Two views of the same question.
%
%   left  -- cumulative weighted variance of eta explained by modes 1..m.
%            Where this flattens, extra modes are buying you nothing.
%   right -- rms residual at each depth after subtracting modes 1..m. If the
%            residual is still structured (not flat, not small) at m = nbc,
%            there is signal you are not capturing.
%
% This is Alford's "know when you look at it", made into something you can put
% on a slide. Adjust cfg.nbc and rerun from cell 5.

ve = median(fit.varexp, 2, 'omitnan');
figure('Position',[80 80 900 400],'Color','w')
subplot(1,2,1)
  plot(1:cfg.nbc, ve, 'o-', 'LineWidth',1.4); grid on
  xlabel('modes retained'); ylabel('median cumulative variance explained')
  ylim([0 1]); xticks(1:cfg.nbc); title('where does this flatten?')
subplot(1,2,2)
  hold on; box on
  for m = 1:cfg.nbc
      r = eta.eta - fit.A(:,1:m)*ahat_raw(1:m,:);
      plot(std(r, 0, 2, 'omitnan'), F2_h.bkgd.zref, 'o-', 'LineWidth',1.2)
  end
  plot(std(eta.eta,0,2,'omitnan'), F2_h.bkgd.zref, 'k--', 'LineWidth',1.2)
  xlabel('rms residual \eta [m]'); ylabel('z [m]')
  set(gca, 'YDir', 'reverse');

  legend([arrayfun(@(m)sprintf('1-%d',m),1:cfg.nbc,'uni',0) {'no fit'}], ...
         'Location','northeast')

fprintf('median cumulative variance explained: %s\n', mat2str(round(ve',3)));
%% ==== % variance explained vs number of modes ======================
nmax = 10;
Mx   = build_modes(F2_h.bkgd.zfull, F2_h.bkgd.N2full, cfg.lat, nmax, cfg.rho0);
Z    = interp1(Mx.z, Mx.Zeta, F2_h.bkgd.zref);          % modes at the fit levels
w    = fit.w(:);                                   % same weights as fit_eta_modes
sets = {eta.eta, eta_M2};  names = {'raw \eta (all frequencies)', 'M_2-band \eta'};
pct  = nan(nmax, 2);
for s = 1:2
    X = sets{s};  num = zeros(nmax,1);  den = 0;
    for it = 1:size(X,2)
        e = X(:,it);  ok = isfinite(e);
        if nnz(ok) < nmax + 2, continue; end        % same hours for every m
        sw  = sqrt(w(ok));
        den = den + sum(w(ok).*e(ok).^2);
        for m = 1:nmax
            a = (sw.*Z(ok,1:m)) \ (sw.*e(ok));       % true m-mode refit
            num(m) = num(m) + sum(w(ok).*(e(ok) - Z(ok,1:m)*a).^2);
        end
    end
    pct(:,s) = (1 - num/den);
end

figure('Color','w')
plot(1:nmax, pct(:,1), 'o-', 'LineWidth',1.4); hold on
plot(1:nmax, pct(:,2), 'o-', 'LineWidth',1.4)
grid on; ylim([0 1]); xlim([0.5 nmax+0.5]); xticks(1:nmax)
xlabel('number of modes fit'); ylabel('median cumulative \eta variance explained')
legend(names, 'Location','southeast')
title('where does this flatten enough?')
fprintf('%% explained, raw: %s\n', mat2str(round(pct(:,1)',1)));
fprintf('%% explained, M2:  %s\n', mat2str(round(pct(:,2)',1)));

%% ==== 7. p', two routes, cross-checked ==============================
% ROUTE A (modal). M.P already carries the sign and the zero-depth-mean
% condition (Nash eq. 4), so this is one line.
pprime_A = F2_h.modes.P * ahat_M2;

% ROUTE B (direct hydrostatic). Rebuild full-depth eta from the modes, convert to
% rho', integrate, demean. Algebraically identical to A -- this checks the
% integration and the sign, not the physics.
eta_full = F2_h.modes.Zeta * ahat_M2;
rho_p    = (F2_h.modes.rho0/F2_h.modes.g) * (F2_h.modes.N2 .* eta_full);
pprime_B = F2_h.modes.g * cumtrapz(F2_h.modes.z, rho_p);
pprime_B = pprime_B - trapz(F2_h.modes.z, pprime_B)/F2_h.modes.H;

d = pprime_A - pprime_B;
fprintf('route A vs B: max |diff| %.3g Pa on %.3g Pa signal (%.1e relative)\n', ...
        max(abs(d(:))), max(abs(pprime_A(:))), max(abs(d(:)))/max(abs(pprime_A(:))));
fprintf('depth-mean of p'' (Nash eq. 4, want 0): max %.3g Pa\n', ...
        max(abs(trapz(F2_h.modes.z, pprime_A)/F2_h.modes.H)));

% Sign check, point by point: dp'/dd = rho0*N2*eta (d positive down, eta positive up).
% Isopycnals raised -> denser water at each depth -> p' increases with depth.
% diff gives dp'/dd at the midpoints between grid levels, so rho0*N2*eta is averaged
% onto the same midpoints. Since P came from cumtrapz, the two should agree to
% round-off, not just in sign.
eta_full = F2_h.modes.Zeta * ahat_M2;
zm   = (F2_h.modes.z(1:end-1) + F2_h.modes.z(2:end))/2;         % midpoints [m]
lhs  = diff(pprime_A, 1, 1) ./ diff(F2_h.modes.z);               % dp'/dd [Pa/m]
f    = F2_h.modes.rho0 * F2_h.modes.N2 .* eta_full;              % rho0*N2*eta [Pa/m]
rhs  = (f(1:end-1,:) + f(2:end,:))/2;

ok   = isfinite(lhs) & isfinite(rhs) & abs(rhs) > 1e-6*max(abs(rhs(:)));
fprintf('sign check: dp''/dd matches sign of rho0*N2*eta at %.2f%% of (depth, time) points\n', ...
        100*mean(sign(lhs(ok)) == sign(rhs(ok))));
fprintf('            max |dp''/dd - rho0*N2*eta| = %.2e of max |rho0*N2*eta|\n', ...
        max(abs(lhs(ok) - rhs(ok)))/max(abs(rhs(ok))));

% one timestep, by eye: the timestep with the largest M2 pressure signal
it=30;
figure
plot(eta_full(:,it), F2_h.modes.z, 'LineWidth',1.4); hold on
plot(pprime_A(:,it), F2_h.modes.z, 'LineWidth',1.4)
xline(0,'k:'); set(gca,'YDir','reverse'); box on
xlabel('\eta (m),  p'' (dbar)')
legend('\eta_{M2}','p''','Location','southeast','Box','off')
title('p'' slopes positive where \eta > 0','FontWeight','normal')


%% ==== 8. the movie ==================================================
fitM = fit;  fitM.varexp = nan(cfg.nbc, numel(F2_h.time));
for it = 1:numel(F2_h.time)
    e = eta_M2(:,it); ok = isfinite(e);
    if nnz(ok) < cfg.fit_minobs || ~all(isfinite(ahat_M2(:,it))), continue; end
    v0 = sum(fit.w(ok).*e(ok).^2);
    for m = 1:cfg.nbc
        r = e(ok) - fit.A(ok,1:m)*ahat_M2(1:m,it);
        fitM.varexp(m,it) = 1 - sum(fit.w(ok).*r.^2)/v0;
    end
end

mopt = struct('outfile','modes_F2_M2_again.mp4', 'fps',8, 'pplot',true);
i0 = find(all(isfinite(ahat_M2),1), 1) + 200*24;    % clear of the filter edges
mopt.frames = i0:2:(i0 + 10*24);                    % 10 days, every 2 h
mode_movie(eta_M2, F2_h.bkgd.zref, F2_h.modes, ahat_M2, fitM, F2_h.time, mopt);


%% ===================================================================
%  F2 MOORING -- CONTEXT THAT IS NOT RECOVERABLE FROM THE DATA ALONE
%  8 30.164'N 134 35.508'E, H = 1515 m, 43 sensors
%  5 Jun 2016 -- 13 Apr 2017 (7512 hourly, from 10-min)
%
%  ------------------------------------------------------------------
%  UNRESOLVED: 22 m DEPTH OFFSET vs THE SPRAY GLIDERS
%  ------------------------------------------------------------------
%  The mooring's thermocline sits ~22 m SHALLOWER than the gliders' over
%  85-250 m; the offset collapses to zero below 300 m. At a fixed depth the
%  mooring reads COLDER. Ruled out so far:
%    - not the wire extrapolation: SBE37/SBE39 sensors with their OWN pressure
%      show the same +21 to +24 m as the RBRs whose depths are interpolated
%    - not time sampling: survives restricting to the glider window (30 Jun
%      onward) and to a single Aug fortnight with a glider inside 15 km
%    - not a wake gradient: glider transects 40 km west and 40 km east agree
%      with EACH OTHER and both disagree with F2 between them
%    - not a water-mass difference: T-S curves coincide
%  Open question for whoever knows the Spray format: is the netCDF 'depth'
%  coordinate (int 10:10:1000) a pressure-derived depth, or binned pressure
%  treated as depth? Emailed Alford, DATE.
%  CONSEQUENCE: do NOT splice glider sigma(z) onto the mooring to fill the
%  upper cap. Anything using glider DEPTH is suspect. Anything using only the
%  glider T-S RELATION is fine -- that is depth-independent and validated.
%
%  ------------------------------------------------------------------
%  SENSOR COVERAGE (checked 9 Sep 2026)
%  ------------------------------------------------------------------
%  row 1 ( 87 m): depth 100%, T 100%, S   0%  <- SBE37, no salinity at all
%  row 2 ( 91 m): depth  41%, T  41%, S   4%  <- dies hour 3049 (~10 Oct 2016)
%  row 3 (107 m): depth 100%, T 100%, S 100%  <- first usable density
%  readme says only: "S: interpolated from glider data to 10 minute
%  intervals". S has 8985 unique values at EVERY depth 107-1190 m, so it was
%  generated by one procedure, not by 43 conductivity cells (there are 4).
%  Whatever that procedure was, it stopped at 107 m.
%  TRAP: zref is built from mean(Dh,2), i.e. DEPTH coverage. Rows with no
%  salinity still get a level, and sig_t there is then invented by
%  fillmissing, not measured. Build zref from sigma coverage instead.
%
%  ------------------------------------------------------------------
%  KNOCKDOWN: REAL, BUT NOT THE EXPLANATION FOR ANYTHING
%  ------------------------------------------------------------------
%  range(depth,2): min 4.1 m, median 56.0 m, max 62.9 m. Episodic (Jul-Aug,
%  Dec-Jan), not continuous. Median profile ~uniform to 260 m then ramps to
%  zero at the anchor.
%  TRAP: "blowdown %.0f m" computed as prctile(ztop,0) - min(ztop) is
%  identically zero by construction. Use range(depth,2).
%  Isopycnal-following eta is immune to knockdown by design -- sigma is
%  interpolated at each sensor's INSTANTANEOUS depth.
%
%  ------------------------------------------------------------------
%  WHY SHALLOW LEVELS GO MISSING: OUTCROPPING, NOT KNOCKDOWN
%  ------------------------------------------------------------------
%  When the 87 m target is absent, sigma at the shallowest sensor is 24.706
%  against a target of 24.247 -- the surface is ABOVE the array 99% of those
%  times. Top sensor depth 104 m present vs 105 m missing: the mooring has
%  not moved, the water changed.
%  So shallow gaps are SELF-CENSORING (seen only when deep), mirroring the
%  deep gaps (seen only when shallow). Both bias the mean one-sidedly.
%  Coverage alone cannot distinguish censoring from harmless intermittency.
%  Use zbar_off = mean(ziso) - zref, NORMALISED by rms eta, since rms eta is
%  ~8 m at the top and ~50 m at depth: |bias| > ~0.3 sigma means the level is
%  measuring half an oscillation. Absolute thresholds fail at one end.
%  Rule of thumb: a trackable isopycnal needs mean depth > (shallowest sigma
%  sensor) + 2*rms eta, i.e. ~120 m here.
%  Currently <50% resolved: 87, 107, 1496 m.
%
%  ------------------------------------------------------------------
%  LARGE DEEP DISPLACEMENTS ARE PHYSICS
%  ------------------------------------------------------------------
%  rms eta reaches 57.1 m at 1289 m. SNR is 21-55 at EVERY level, so this is
%  not a resolution limit. rms eta scales as N^-1 to within 3% and APE
%  density is ~constant with depth. Do not trim deep levels for being large;
%  trim them for zbar_off (reaches -66 m below 1400 m).
%
%  ------------------------------------------------------------------
%  THE TIDE IS MODE-2 DOMINANT IN DISPLACEMENT
%  ------------------------------------------------------------------
%  rms M2 amplitudes [7.4 19.0 9.1 12.5] m. NOT continuum leaking through a
%  wide filter: peak-above-continuum is 65-78% of band variance for all four
%  modes. Mode 2 genuinely dominates eta. Consistent with a nearby generation
%  site -- 1515 m of water means this is a slope, not the abyss.
%  Mode 1 still leads p' because p' is a vertical integral of xi and so is
%  redder by k_z^2 (Nash 2005). Modes 1-2 carry ~79% of p' but only ~54% of
%  M2 displacement variance.
%  Variance explained against RAW eta (0.665 at 4 modes) is the WRONG number
%  -- against M2-band eta it is ~0.86. Raw eta includes the broadband
%  continuum, which four modes were never going to represent.
%  Diurnal matters here: K1/O1 peak rivals M2, and the K1 critical latitude
%  (30 deg) is far away, so diurnal internal tides propagate freely. Running
%  the pipeline at K1 costs one filter call.
%
%  ------------------------------------------------------------------
%  MIXED LAYER: THE LAST OPEN ASSUMPTION
%  ------------------------------------------------------------------
%  Nothing measures density above the shallowest sigma sensor. c1 spans
%  1.80-2.00 m/s across treatments of that cap; current run gives
%  [1.878 1.017 0.614 0.443] m/s. Mode 1's rms p' spikes to ~150 Pa at the
%  surface, INSIDE the assumed layer, so this is not a harmless boundary
%  condition -- it shapes the quantity the flux is built on.
%  A neutral cap is conservative in asserting nothing, but it is not neutral
%  in effect: c ~ (1/pi)*int N dz, so removing stratification biases c1 low.
%  Report the spread rather than hiding it.
%
%  ------------------------------------------------------------------
%  sw_vmodes GOTCHAS (verified numerically)
%  ------------------------------------------------------------------
%  1. Column 1 is the BAROTROPIC mode (PVel(1) = sqrt(gH), Hori(:,1) flat to
%     0.5%). Ask for nbc+1, take columns 2:nbc+1.
%  2. Vert and Hori are each separately renormalised to max = 1 on the last
%     two lines, which destroys the amplitude link between them. Rebuild
%     pressure modes by integrating N^2*Vert directly.
%  3. Vert is NOT orthogonal under int N^2 Zi Zj once column 1 is included
%     (cols 1,2 correlate -0.88). Baroclinic columns alone are orthogonal.
%  4. It imposes a FREE SURFACE, not a rigid lid. Baroclinic Vert(1,:) comes
%     back at ~0.003, so Alford's eq. A2 rigid-lid statement is fine to
%     quote; the difference only matters for the mode being discarded.
%
%  ------------------------------------------------------------------
%  CHECKS THAT SHOULD PASS EVERY RUN
%  ------------------------------------------------------------------
%  route A vs B p':      ~1e-15 relative           (currently 2.9e-15)
%  depth-mean of p':     ~0 Pa, Nash eq. 4         (currently 4e-14)
%  sign: eta > 0 pushes isopycnals up, adding mass above, so interior p'
%        goes POSITIVE  
%  max |time-mean eta|:  ~0 m                      (currently 0.31)
%  fit_first vs filter_first: a few % if healthy   (currently 8-14%)
%  Thorpe sorting fires on 29% of profiles, median 2 levels -- high; this is
%  an assertion that overturns in a 43-sensor hourly profile are artifacts.
%% ===================================================================

% Same baroclinic set as the displacement fit, plus mode 0. Passing nmodes
% explicitly is the place to narrow the velocity fit if the ADCP count says to
% (see the degrees-of-freedom note in fit_uv_modes); cfg.nbc is the default.
F2_h.uv_modes = fit_uv_modes(F2_h, F2_h.modes, cfg, struct('nmodes', cfg.nbc));

nrow = size(F2_h.uv_modes.uhat_M2, 1);   % 1 + however many modes were actually fit
c    = orderedcolors("gem");

figure('Position',[80 80 1000 650],'Color','w')
for r = 1:nrow
    subplot(nrow, 1, r)
    plot(F2_h.time, F2_h.uv_modes.uhat_M2(r,:), 'Color', c(1,:)); hold on
    datetick('x','mmm'); box on
    if r == 1
        ylabel('u m0 (BT) [m/s]')
    else
        ylabel(sprintf('u m%d [m/s]', r-1))
    end
end
linkaxes(findall(gcf,'Type','axes'),'y')

figure('Position',[80 80 1000 650],'Color','w')
for r = 1:nrow
    subplot(nrow, 1, r)
    plot(F2_h.time, F2_h.uv_modes.vhat_M2(r,:), 'Color', c(2,:))
    datetick('x','mmm'); box on
    if r == 1
        ylabel('v m0 (BT) (m/s)')
    else
        ylabel(sprintf('v m%d (m/s)', r-1))
    end
end
linkaxes(findall(gcf,'Type','axes'),'y')

figure('Position',[80 80 1000 650],'Color','w')
for r = 1:nrow
    subplot(nrow, 1, r)
    plot(F2_h.time, F2_h.uv_modes.uhat_M2(r,:), 'Color', c(1,:)); hold on
    plot(F2_h.time, F2_h.uv_modes.vhat_M2(r,:), 'Color', c(2,:))
    datetick('x','mmm'); box on
    if r == 1
        ylabel('mode 0 (BT) (m/s)')
        legend('u','v','Box','off','Orientation','horizontal')
    else
        ylabel(sprintf('mode %d [m/s]', r-1))
    end
end
linkaxes(findall(gcf,'Type','axes'),'y')


%% tidal flux calculation
tidal_flux = compute_tidal_flux(ahat_M2, F2_h.uv_modes.uhat_M2, F2_h.uv_modes.vhat_M2, ...
                       F2_h.modes, cfg, struct('periods', 2));

flux_movie(tidal_flux, F2_h.time, 8.5027, 134.5918, ...
           'Palau_4_100m_100m_139merge.nc', struct('outfile','flux_F2.mp4'));


%% tidal ellipse to check...

E = plot_tidal_ellipse(F2_h, cfg);

%%
R = synth_check(F2, F2_h, eta, ahat_M2, tidal_flux, cfg);