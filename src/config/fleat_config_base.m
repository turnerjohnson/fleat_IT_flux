function cfg = fleat_config_base(name)
% FLEAT_CONFIG_BASE  Shared defaults for every FLEAT mooring config.
%
%   cfg = fleat_config_base('F3')
%
% WHAT THIS FILE IS FOR
% f2_config and f4_config were ~180 lines each and about 90% the same text: the
% same defaults with the same paragraph of reasoning beside each one. Five
% moorings would have been five copies, and a default that changed would have
% had to change in five places. Everything that is the SAME for every mooring
% lives here once, together with the reasoning. Each fN_config then sets only
% what genuinely differs at that site, and the diff between two moorings is
% readable in a screen.
%
% The rule for deciding where something goes:
%   here          a default, or a statement about the method
%   in fN_config  an assumption or a measurement about THAT mooring
%
% So mld, the record window, the T/S gates and anything marked unverified
% belong to the mooring. nbc, the bandwidth, the smoothing window and the
% filtering order belong to the method, and changing one here changes it for
% the whole array, which is what you want when you are comparing moorings.
%
% author TEJ, with help from Claude
% 15 October 2026

%% ---- site ------------------------------------------------------------
% THE POSITION TABLE FOR THE WHOLE ARRAY, in degrees and decimal minutes off
% the deployment mooring table, converted by dm() at the bottom rather than
% retyped in decimal degrees.
%
% This used to be in two places: POS in fleat_moorings for anything looping
% over the array, and a private copy in f2_config and f4_config so a
% single-mooring run needed nothing but its config. Both of those are now
% served from here, so there is one set of numbers. fleat_moorings still has
% its own POS; when you next touch that file, have it read these instead.
%
% lat is NOT in the data. See the note in fleat_moorings: info.lat in the
% deployment file holds LONGITUDE at every mooring, and there is no latitude
% anywhere in the file, the fullfield products or the instrument lists. These
% numbers are the only record of it.
POS = struct( ...
    'F1', [dm(8,40.838)  dm(134,24.087)], ...
    'F2', [dm(8,30.164)  dm(134,35.508)], ...
    'F3', [dm(8,32.723)  dm(134,37.414)], ...
    'F4', [dm(8,29.862)  dm(134,38.723)], ...
    'F6', [dm(8,28.962)  dm(134,38.575)]);

% Recorded bottom depth [m].
HREC = struct('F1', 3390, 'F2', 1515, 'F3', 1666, 'F4', 880, 'F6', 434);

name = upper(string(name));
name = char(name);
if ~isfield(POS, name)
    error('fleat_config_base:unknown', ...
        ['"%s" is not in the position table. Known moorings: %s. Add it to ' ...
         'POS and HREC in this file if a sixth deployment turns up.'], ...
        name, strjoin(fieldnames(POS)', ', '));
end

cfg.mooring = name;
cfg.lat     = POS.(name)(1);
cfg.lon     = POS.(name)(2);
cfg.H       = HREC.(name);      % water depth [m]
cfg.rho0    = 1025;             % reference density [kg m^-3]

% making sure that we have the needed inputs:
cfg.unverified = {};
if ~isfinite(cfg.H), cfg.unverified = [cfg.unverified {'H'}]; end

%% ---- data location ---------------------------------------------------
% T, S, z from Anna Simpson; u, v from Conrad Luecke.
cfg.datadir  = '/Users/tej/Documents/research_fleat/data/simpson/';
cfg.datafile = fullfile(cfg.datadir, sprintf('%s_fullfield.mat', name));
cfg.datavar  = 'fullfield';     % the variable inside the .mat

% Bathymetry for the map plots (flux_movie, plot_tidal_ellipses_latlon). 
cfg.bathy    = '/Users/tej/Documents/research_fleat/data/Palau_4_100m_100m_139merge.nc';

%% ---- time ------------------------------------------------------------
cfg.avg_k   = 6;                % 10-min samples per hourly average
cfg.dt_hr   = 1;                % resulting timestep [h]

%% ---- record window ---------------------------------------------------
% [] on either end keeps everything, which is what trim_mooring_record tests
% for. A mooring that loses an instrument mid-record overrides these; F4 is the
% example, and the reason lives in f4_config next to the date.
%
% This trims BEFORE the hourly average, so the background mean, N2, the modes
% and every fit all see the same window.
cfg.t_start = [];               % datenum, or [] for the start of the record
cfg.t_end   = [];               % datenum, or [] for the end of the record

%% ---- T/S quality control ---------------------------------------------
% Applied to the raw samples, before the hourly average (see qc_ts). [] on
% either disables that test.
%
% WHY ANY MOORING NEEDS THIS: sw_dens0 evaluates S.^1.5, so one negative
% salinity returns a COMPLEX sigma, and sort, unique and the relational
% operators all accept complex input without complaining. Nothing fails until
% an interp1 several functions downstream. F4's gridded product carries
% salinities from -114 to +489 psu at its top five sensors; F2's never leave
% [34.33 35.07].
%
% The default is the guard rather than the trust. [0 33] and [33 36] are F4's
% values; at F2 they remove nothing at all, so switching them on everywhere
% costs that mooring nothing and would catch a regridding that went wrong at
% one of the three moorings nobody has looked at yet.
%
% A hard range is deliberate. A robust per-sensor test (median +/- n*MAD)
% rejects ten times as much at the shallow sensors, where the spread is real
% upper-ocean variability rather than instrument failure.
cfg.T_valid = [0 33];           % [degC] [lo hi], or [] to disable
cfg.S_valid = [33 36];          % [psu]  [lo hi], or [] to disable

%% ---- frequency bands -------------------------------------------------
% f_cpd is derived from cfg.lat and so is set here, after the position, where
% it cannot go stale. If you ever override cfg.lat in an fN_config, recompute
% this on the next line or it will describe the wrong site.
cfg.f_cpd   = 2*sind(cfg.lat);  % 0.296 cpd across the array. Inertial 3.4 d.
cfg.M2_cpd  = 1.9323;
cfg.K1_cpd  = 1.0027;
cfg.bw_c    = 1.25;             % Alford (2003) A2 bandwidth parameter


% Diurnal matters here too: the K1/O1 peak rivals M2 at F2 and the K1 critical
% latitude (30 deg) is far away, so diurnal internal tides propagate freely.
% Running the array at K1 costs one field: set cfg.band_cpd below.

% WHICH BAND THE PIPELINE ACTUALLY RUNS. Named separately from M2_cpd so that a
% K1 run is one line in the driver (cfg.band_cpd = cfg.K1_cpd) rather than an
% edit in several places.
%
% HOW IT REACHES THE REST OF THE CODE: bandpass_alford is called with
% cfg.M2_cpd inside fit_uv_modes, and compute_tidal_flux reads cfg.M2_cpd to
% set its phase-average window, so cfg.M2_cpd is the field those functions
% actually consume. fleat_IT_flux copies band_cpd onto M2_cpd on its first
% line and says so, which gets a K1 run without editing either of them. The
% *_M2 variable names downstream are therefore historical: after that copy they
% mean "the band being run".
cfg.band_cpd  = cfg.M2_cpd;
cfg.band_name = 'M2';

%% ---- order of filtering and fitting: SETTLED, no longer a config knob -
% There used to be a cfg.band_order here, 'fit_first' | 'filter_first'. It was
% removed on 10 Nov 2026 and the pipeline now always fits first.
%
% Alford does filter_first: Appendix A3 defines eta from "the band-passed
% temperature at depth z_j", and A4 fits modes to "the filtered time series...
% at discrete depths z_j". Nash states the two orders are equivalent, noting
% that harmonic analysis may equivalently be performed on time-series of
% vertical mode amplitudes.
%
% They are exactly equivalent only if the fit is the same linear operator at
% every timestep, which fails slightly when a level drops out and the design
% matrix loses a row. fit_first wins for a practical reason: filtfilt is
% recursive and cannot take a NaN, so one order or the other has to fill gaps
% first -- filter_first on ~43 gappy sensor levels, fit_first on nbc smooth
% amplitude series. That is much less invented data, at every mooring in this
% array, so there was nothing left for the knob to decide.
%
% Alford's order survives as a DIAGNOSTIC: opt.check_order in fleat_IT_flux
% refits filter_first and prints the relative difference by mode, 8-14% at F2.
% Run it once per mooring; if the amplitudes agree the choice stops mattering
% and you can say so in a sentence. The reasoning lives at the top of cell 5 in
% fleat_IT_flux, which also errors if an old config still sets band_order, so a
% stale filter_first cannot run silently.

%% ---- background density profile --------------------------------------
cfg.zref_mode = 'instrument';   % 'instrument' | 'uniform'
cfg.zref_dz   = 10;             % [m] grid spacing if zref_mode = 'uniform'
cfg.zref_merge= 8;              % [m] merge sensors closer than this ('instrument')
cfg.sig_cover_min = 0.2;        % sensor must report sigma this fraction of hours
                                % to contribute a background level

% MIXED LAYER: the last open assumption, and the one field every new config has
% to think about. Nothing measures density above the shallowest sigma sensor, so
% that cap is asserted, not measured. Setting mld at or below the shallowest
% zref level makes the whole unmeasured cap neutral, which is conservative in
% asserting nothing -- but it is NOT neutral in effect, because
% c ~ (1/pi)*int N dz, so removing stratification biases c1 low. At F2, c1 spans
% 1.80-2.00 m/s across treatments of the cap. Report the spread.
%
% NaN here forces each mooring to state it. Set it from that mooring's
% shallowest sigma-bearing sensor, which f_background_v2 prints.
% fleat_IT_flux refuses to run on a NaN mld rather than letting it reach
% f_background_v2, so this is a required field, not an unverified one.
cfg.mld      = NaN;             % [m]


cfg.mld_dsigma   = 0.03;        % [kg m^-3] density-threshold criterion
cfg.mld_pea_crit = 10;          % [J m^-2] PEA criterion, or [] (see above)
cfg.mld_zref     = [];          % [m] fixed reference depth, or [] for hourly
cfg.mld_zmax     = 400;         % [m] how far below the reference to search
cfg.mld_dz       = 1;           % [m] grid phi(h) is evaluated on
cfg.mld_minobs   = 4;           % density samples needed to try an hour

% The index to use when the mixed layer base is NOT resolved, which here is
% most hours: phi over a fixed thickness below the reference. It needs no base
% and so is defined in every hour with a profile, which makes it the honest
% thing to plot against time and to compare across the array. Reichl's phi
% with h prescribed rather than solved for, in other words.
cfg.mld_phi_span = 100;         % [m]



cfg.N2_floor = 1e-7;            % [s^-2] sw_vmodes rejects N <= 0

%% ---- which N2 --------------------------------------------------------
% TWO INDEPENDENT ESTIMATES OF THE SAME PROFILE. f_background_v2 builds both
% every run and this only decides which one lands in bkgd.N2full and so reaches
% build_modes.
%
%   'sigma'    difference the record-mean S and T at the reference levels with
%              sw_bfrq. Ours end to end: we choose the levels, the rejection,
%              the smoothing and the mixed-layer cap, and all of it is visible
%              in f_background_v2.
%   'gridded'  the time median of the N2 field Simpson supplies in F.N2, at the
%              mean sensor depths. Owes nothing to any of those choices, and
%              reaches sensors that have no salinity and so no reference level
%              of their own.
%
% Neither is the truth. Look at the middle panel of plot_background, then flip
% it with set_N2_source on bkgd -- both profiles are already there, so the
% switch costs nothing and does not rebuild the background. If the eigenspeeds
% move a lot between the two, that spread is a result worth reporting, not a bug
% to be tuned away.
cfg.N2_source = 'sigma';        % 'sigma' | 'gridded'

% Gridded branch only. Above the shallowest gridded level there is no estimate,
% so N2 takes cfg.N2_ml, the same neutral cap the sigma branch uses. Below the
% deepest, a least-squares straight line over the deepest N2_grid_fit_m metres
% is continued to cfg.H and floored at cfg.N2_floor.
% A TWO-POINT slope would be identically zero: the readme says dTdz is a
% "forward average between successive instruments", so F.N2 repeats its deepest
% value, and at F2 and F4 alike the last two levels are equal to all printed
% digits. Hence a span rather than the last pair.
cfg.N2_grid_fit_m = 100;        % [m] depth span of the bottom linear fit

%% ---- N2 smoothing ----------------------------------------------------
% Smoothed in log space because N2 spans three decades; linear smoothing would
% flatten the thermocline peak. Matthew's "smooth nsq" note.
%
% N2_smooth_m is f_background_v2's window and it is in METRES, not points,
% because the reference levels run 10 to 148 m apart and a fixed point count
% would be a 50 m filter in the thermocline and a 500 m filter through the
% middle of the column. On a uniform grid it reduces exactly to Shaun's
% nanrunmean(Nsq, nzsmooth).
cfg.N2_smooth_m    = 50;        % [m] vertical running-mean window
cfg.N2_ml          = cfg.N2_floor;  % [s^-2] N2 in the unmeasured cap

% Read only by the older f2_background. Left here so that function still runs.
cfg.N2_smooth_span = 9;
cfg.N2_smooth_meth = 'loess';
cfg.zfull_dz       = 2;         % [m] mode solution grid

%% ---- displacement ----------------------------------------------------
cfg.eta_ref_days = 25;          % running-mean window defining the reference depth
cfg.min_levels   = 6;           % min good sensors to attempt an isopycnal profile
cfg.min_coverage = 0;           % drop zref levels resolved less often than this

%% ---- modes -----------------------------------------------------------
% One mode count for the whole pipeline AND the whole array. Displacement and
% velocity are fit to the SAME baroclinic set: the flux sum collapses to
% min(n_pressure, n_velocity) anyway (compute_tidal_flux), so a separate count
% for each side bought nothing. To fit velocity to FEWER modes than eta -- which
% the degrees-of-freedom note in fit_uv_modes argues for, given three ADCPs --
% pass opt.uv_nmodes to fleat_IT_flux rather than adding a second config field.
%
% Keeping one nbc across the array is a choice worth being deliberate about. The
% moorings span 434 to 1515 m of water, so mode 6 at F6 is a much shorter
% vertical wavelength than mode 6 at F2 and is resolved by far fewer sensors.
% The honest check is per-mooring: look at opt.nmode_sweep and at the varexp
% diagnostic, and if a shallow mooring cannot support six, say so and drop ITS
% nbc in its own config rather than lowering the array.
cfg.nbc          = 6;             % baroclinic modes (velocity adds mode 0 in the fit)
cfg.fit_weight   = 'wkb';         % 'wkb' | 'none'
cfg.fit_minobs   = 8;             % min good levels to attempt a fit at a timestep
cfg.fit_maxcond  = 100;           % reject a timestep whose weighted design matrix
                                  % is worse conditioned than this (Inf disables).
% minobs counts levels; maxcond asks whether they span the column. Without it,
% 39 of 3936 hours at F4 -- hours when the thermocline drops out and only deep
% levels survive -- return amplitudes of 20000 to 39000 m and set the scale of
% every amplitude plot. See the note in fit_eta_modes.

%% ---- flux ------------------------------------------------------------
% Phase-average window for <u'p'>, in tidal periods. A boxcar over an integer
% number of periods removes the 2w term exactly, but M2 is 12.4206 h and the
% sampling is hourly, so 2 periods is 24.84 h rounded to 25 and leaves ~0.6%
% leakage. More periods reduces it at the cost of temporal resolution; see
% compute_tidal_flux, and opt.flux_lowpass to sidestep the rounding entirely.
cfg.flux_periods = 2;
end

%% ======================================================================
function dd = dm(d, m)
% degrees and decimal minutes -> decimal degrees. Northern and eastern
% hemisphere only, which is where this array is.
dd = d + m/60;
end
