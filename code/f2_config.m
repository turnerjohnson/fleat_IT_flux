function cfg = f2_config()
% F2_CONFIG  Every choice in the F2 modal-decomposition pipeline, in one place.

% Several of these are assumptions rather than measurements (mld, N2 smoothing,
% zref_mode, nbc). Having them in one struct makes it a two-line job to rerun the
% pipeline with a different value and show how much the answer moved.

% ---- site -----------------------------------------------------------
cfg.mooring = 'F2';
cfg.lat     = 8.5;              % 08 30.164' N
cfg.H       = 1515;             % water depth [m]
cfg.rho0    = 1027;             % reference density [kg m^-3]

% ---- time -----------------------------------------------------------
cfg.avg_k   = 6;                % 10-min samples per hourly average
cfg.dt_hr   = 1;                % resulting timestep [h]

% ---- frequency bands ------------------------------------------------
cfg.f_cpd   = 2*sind(cfg.lat);  % 0.2957 cpd at F2. Inertial period 3.4 days.
cfg.M2_cpd  = 1.9323;
cfg.K1_cpd  = 1.0027;
cfg.bw_c    = 1.25;             % Alford (2003) A2 bandwidth parameter

% At 8.5 N the near-inertial and semidiurnal bands are 1.64 cpd apart. Alford's
% warning that they become inseparable poleward of 50 deg does not apply, and the
% c = 1.25 filter isolates M2 cleanly. A real advantage of F2 over the
% midlatitude moorings in the Internal Swell paper.

% ---- order of filtering and fitting ---------------------------------
% 'fit_first'    fit raw eta, then bandpass the modal amplitude series
% 'filter_first' bandpass eta at the 43 levels, then fit

% Alford does filter_first: Appendix A3 defines eta from "the band-passed
% temperature at depth z_j", and A4 fits modes to "the filtered time series...
% at discrete depths z_j". Nash states the two orders are equivalent, noting
% that harmonic analysis may equivalently be performed on time-series of
% vertical mode amplitudes.

% They are exactly equivalent only if the fit is the same linear operator at
% every timestep, which fails slightly when a level drops out and the design
% matrix loses a row. fit_first is the default here for a practical reason:
% filter_first needs gap-filling on 43 series before filtfilt, fit_first needs
% it on nbc. That is much less invented data. Run both once and compare -- if
% the amplitudes agree, the choice stops mattering and you can say so.
cfg.band_order = 'fit_first';

% ---- background density profile -------------------------------------
cfg.zref_mode = 'instrument';   % 'instrument' | 'uniform'
cfg.zref_dz   = 10;             % [m] grid spacing if zref_mode = 'uniform'
cfg.zref_merge= 8;              % [m] merge sensors closer than this ('instrument')
cfg.sig_cover_min = 0.2;        % sensor must report sigma this fraction of hours
                                % to contribute a background level                                                % prev 85

% Mixed layer above the shallowest sensor (Alford's recommendation). The top
% T sensor is at 78 m and nothing is above it, so the upper 78 m of the mean
% profile is asserted, not measured. Setting mld at or below the shallowest zref
% level makes the whole unmeasured cap neutral, which is the conservative
% choice: no artificial stratification where there is no data.

% If the FLEAT gliders profiled the upper 78 m, use them here instead and this
% assumption goes away.
cfg.mld      = 117;              % [m]
cfg.N2_floor = 1e-7;            % [s^-2] sw_vmodes rejects N <= 0

% ---- N2 smoothing ---------------------------------------------------
% Smoothed in log space because N2 spans three decades; linear smoothing would
% flatten the thermocline peak. Matthew's "smooth nsq" note.
cfg.N2_smooth_span = 9;
cfg.N2_smooth_meth = 'loess';
cfg.zfull_dz       = 2;         % [m] mode solution grid

% ---- displacement ---------------------------------------------------
cfg.eta_ref_days = 25;          % running-mean window defining the reference depth
cfg.min_levels   = 6;           % min good sensors to attempt an isopycnal profile
cfg.min_coverage = 0;         % drop zref levels resolved less often than this

% ---- modes ----------------------------------------------------------
% One mode count for the whole pipeline. Displacement and velocity are fit to
% the SAME baroclinic set: the flux sum collapses to min(n_pressure, n_velocity)
% anyway (compute_tidal_flux), so a separate count for each side bought nothing.
% To fit velocity to FEWER modes than eta -- which the degrees-of-freedom note
% in fit_uv_modes argues for, given three ADCPs -- pass opt.nmodes explicitly at
% that call rather than reintroducing a second config field.
cfg.nbc          = 6;             % baroclinic modes (velocity adds mode 0 in the fit)
cfg.fit_weight   = 'wkb';         % 'wkb' | 'none'
cfg.fit_minobs   = 8;             % min good levels to attempt a fit at a timestep
end