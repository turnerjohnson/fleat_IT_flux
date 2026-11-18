function cfg = f4_config()
% F4_CONFIG  What is true of the F4 mooring, on top of fleat_config_base.
%
%   8 29.862'N 134 38.723'E, H = 880 m, 34 instruments (56 to 848 m)
%   Record to 13 Apr 2017; THIS CONFIG USES 164 of 313 days (see below).
%
% Everything shared with the rest of the array, and the reasoning behind it, is
% in fleat_config_base. This file is only F4.
%
% T, S, z from Anna Simpson; u, v from Conrad Luecke.

cfg = fleat_config_base('F4');

%% ---- record window ---------------------------------------------------
% THE DEEP ADCP FAILS ON 2016-11-18. Sensors 25-33 (640 to 824 m) report u and
% v with full coverage up to 2016-11-17 and 1% of the time after it, so every
% velocity fit past that date extrapolates the bottom two thirds of the column
% from the top 450 m. That is where the amplitudes go wrong.
%
% Two other things change at roughly the same point and argue the same way: the
% array settles about 25 m deeper, and sensor 2 loses temperature from
% 2017-02-04. Cutting at 736650 (2016-11-16) stops two days short of the ADCP
% failure rather than trying to ride it out.
%
% CONSEQUENCE FOR THE ARRAY COMPARISON: F4's flux record is 164 days against
% F2's 313. Any array-mean or any statement about seasonality has to say so,
% and the movie has to tolerate F4 going blank partway through. align_flux_records
% handles that by padding with NaN rather than by shortening everyone to F4.
cfg.t_end   = 736650;           % 2016-11-16. 164 of 313 days.

%% ---- T/S quality control ---------------------------------------------
% F4 IS WHY THE GATES IN THE BASE EXIST. Its gridded product carries salinity
% spikes from -114 to +489 psu at the top five sensors, and a negative S makes
% sw_pden return a COMPLEX sigma that travels silently all the way to the
% interp1 inside f2_eta. Sensors below ~110 m never leave [34.38 35.07], so the
% window only bites at the top and costs 257 of 1.46 million samples (0.018%).

%% ---- mixed layer -----------------------------------------------------
% Shallowest sigma-bearing sensor. Same argument as F2, and the same caveat:
% the cap is asserted, not measured, and it biases c1 low.
cfg.mld = 50;                   % [m]

%% ---- modes -----------------------------------------------------------
% nbc = 6 and fit_maxcond = 100 both come from the base. maxcond is not
% cosmetic here: without it, 39 of 3936 hours at F4 -- hours when the
% thermocline drops out and only deep levels survive -- return amplitudes of
% 20000 to 39000 m and set the scale of every amplitude plot.

end
