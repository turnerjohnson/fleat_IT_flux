function cfg = f1_config()
% F1_CONFIG  What is true of the F1 mooring, on top of fleat_config_base.
%
%   8 40.838'N 134 24.087'E -- the northwestern mooring, furthest from
%   Velasco Reef and furthest from F3/F4/F6.
%
%   *** PLACEHOLDER. NOT YET VERIFIED AGAINST THE DATA. ***
%
% Position comes from the POS table in the base and is as good as F2's: both
% are off the same deployment table. Everything below is what still has to be
% established before an F1 flux means anything. fleat_IT_flux prints
% cfg.unverified in a banner so a run of this config cannot be mistaken for a
% finished one.

cfg = fleat_config_base('F1');

%% ---- TO ESTABLISH, in this order -------------------------------------
%

% lower ADCP failure in Novebver
cfg.t_end   = 736675; 

% 3. SHALLOWEST SENSOR THAT CARRIES SALINITY -> cfg.mld. Not the shallowest
%    sensor: at F2 the top two have no usable S and the first real density is
%    the third one down. f_background_v2 prints the levels it kept.

cfg.mld = 100;                  % [m] REQUIRED -- see step 3. Pipeline errors on NaN.

%
% 5. nbc. With H and the sensor count in hand, run opt.nmode_sweep and see
%    whether six modes are supported. Drop cfg.nbc here if they are not.
%
cfg.nbc = 3;

end
