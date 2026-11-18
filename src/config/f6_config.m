function cfg = f6_config()
% F6_CONFIG  What is true of the F6 mooring, on top of fleat_config_base.

cfg = fleat_config_base('F6');

%% ---- TO ESTABLISH, in this order -------------------------------------
%
% 3. SHALLOWEST SENSOR THAT CARRIES SALINITY -> cfg.mld. Not the shallowest
%    sensor: at F2 the top two have no usable S and the first real density is
%    the third one down. f_background_v2 prints the levels it kept.
cfg.mld = 30;                  % [m] REQUIRED -- see step 3. Pipeline errors on NaN.





%% ---- MODE COUNT: SET BY THE LEVEL COUNT, AND IT IS 7 ----------------
% The gridded product has 8 rows and row 8 (RBR 72141, 419.5 m) carries NO
% salinity and NO velocity -- it hangs on the releases below the deepest ADCP
% bin. f_background_v2 keeps a reference level only where S and T are both
% finite, so F6 has SEVEN density levels and seven velocity levels, not eight.
%
% Against 7 levels, on the "unknowns + 2" rule the fit gates use:
%
%   eta       nbc unknowns    needs nbc+2 <= 7  ->  nbc <= 5
%   velocity  nbc+1 unknowns  needs nbc+3 <= 7  ->  nbc <= 4   <-- binds
%
% Velocity is the limit, not density. The flux sum collapses to
% min(n_pressure, n_velocity) anyway (see compute_tidal_flux and the nbc note
% in fleat_config_base), so asking eta for a 5th mode the velocity side cannot
% supply buys nothing. nbc = 4.
%
% nbc = 4 at nz = 7 is the floor, not a comfortable choice: BOTH fits then need
% all 7 levels in an hour, and f_background_v2 warns about exactly this. Expect
% a large NaN fraction. Drop to 3 if too little of the record survives.
%
% THE COUNT IS NOT THE REAL WORRY. The 7 levels span 137 to 420 m of a 434 m
% column, so the top 137 m -- 32% of the water column, where mode 1 and 2
% displacement structure peaks -- is unsampled entirely. Nothing in the level
% count sees that; cfg.fit_maxcond does. Look at what fraction of hours it
% rejects before believing mode 4 here.
%
% The resolution is recoverable in principle: the deployment sheet lists 19
% SBE 56 thermistors between 131 and 405 m that are not in the gridded product
% (f6_info.txt sec 2). Either they were never gridded or they were lost. That
% question is worth more than any choice of nbc on this file.
cfg.nbc = 4;

% The base value of 8 exceeds F6's level count, so EVERY hour would be rejected
% on both sides and the run would return all-NaN rather than error.
%
% 6, NOT 7. 7 demands that every level report, and the shallowest one does not:
% the eta step reports 137 m resolved less than half the time, which cost 5571
% of 7512 hours and left the fit solving 26% of the record. 6 is the nbc+2 floor
% for nbc = 4, so it asks exactly what the fit needs and no more.
%
% What this buys, and what it costs. The recovered hours are ones where 137 m
% dropped out, so they fit 4 modes across 177-377 m -- 200 m of a 434 m column.
% That is a conditioning question, not a counting one, and cfg.fit_maxcond (100,
% from the base) is what answers it. Expect the "ill-conditioned" count in the
% eta fit line to stop being 0. That is the guard doing its job; if it rejects
% most of the recovered hours, the levels genuinely cannot separate 4 modes and
% nbc = 3 is the answer, not a looser maxcond.
cfg.fit_minobs = 6;

cfg.unverified = [cfg.unverified {'nbc'}]; %'mld','t_end','T_valid','S_valid',
end
