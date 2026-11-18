function cfg = f3_config()
% F3_CONFIG  What is true of the F3 mooring, on top of fleat_config_base.
%
%   8 32.723'N 134 37.414'E -- between F2 and F4, in the eastern cluster.
%
%   *** PLACEHOLDER. NOT YET VERIFIED AGAINST THE DATA. ***
%
% Position comes from the POS table in the base and is as good as F2's: both
% are off the same deployment table. Everything below is what still has to be
% established before an F3 flux means anything. fleat_IT_flux prints
% cfg.unverified in a banner so a run of this config cannot be mistaken for a
% finished one.

cfg = fleat_config_base('F3');

%% ---- TO ESTABLISH, in this order -------------------------------------
%
% 1. BOTTOM DEPTH. Run fleat_moorings(), which reads the recorded depth from
%    the deployment file and checks it against the bathymetry patch; the
%    docstring says F3 lands on its recorded depth exactly. Paste the number
%    into HREC in fleat_config_base, not here, so fleat_moorings and the
%    configs keep agreeing.
%
% 2. THE RECORD AND ITS INSTRUMENT FAILURES. Load the fullfield and look at u,
%    v and T coverage per sensor against time BEFORE anything else. F4's deep
%    ADCP quit on 2016-11-18 and nothing in the pipeline complains about it:
%    the velocity fit happily extrapolates two thirds of the column and returns
%    plausible-looking amplitudes. Set cfg.t_end if F3 does the same.
%      sum(isfinite(F.u),2)/size(F.u,2)   by sensor, and the same for T and S
%
% 3. SHALLOWEST SENSOR THAT CARRIES SALINITY -> cfg.mld. Not the shallowest
%    sensor: at F2 the top two have no usable S and the first real density is
%    the third one down. f_background_v2 prints the levels it kept.
%
% 4. T/S GATES. The base turns them on everywhere ([0 33], [33 36]). Check what
%    qc_ts removes at F3: zero means F3 is clean like F2, a few hundred at the
%    top sensors means it is like F4, and thousands spread down the column means
%    stop and look at the product rather than gating harder.
%
% 5. nbc. With H and the sensor count in hand, run opt.nmode_sweep and see
%    whether six modes are supported. Drop cfg.nbc here if they are not.
%
% Nothing else should need changing. If something does, that is a finding about
% F3 and it belongs in this file next to the number it changed.

cfg.mld = 50;                  % [m] REQUIRED -- see step 3. Pipeline errors on NaN.

cfg.unverified = [cfg.unverified {'mld','t_end','T_valid','S_valid','nbc'}];
end
