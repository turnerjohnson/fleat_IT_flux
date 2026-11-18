function cfg = f2_config()
% F2_CONFIG  What is true of the F2 mooring, on top of fleat_config_base.
%
%   8 30.164'N 134 35.508'E, H = 1515 m, 43 sensors
%   5 Jun 2016 -- 13 Apr 2017 (7512 hourly, from 10-min)
%
% Everything shared with the rest of the array, and the reasoning behind it, is
% in fleat_config_base. This file is only F2: its assumptions, and the
% findings that justify them. The long F2 field notes that used to sit in the
% middle of turner_fleat2_IT_flux are at the bottom.
%
% Position and H come from the POS and HREC tables in the base, so they are not
% retyped here; lat is 8.5027 at full precision, which is the value every
% plotting call used to carry as a literal.

cfg = fleat_config_base('F2');

%% ---- record window ---------------------------------------------------
% F2 has no instrument failure that forces a cut, so it uses its whole record.
% Stated rather than left to the default, because it is a result: the F4 window
% exists for a reason and F2's absence of one is worth the line.
cfg.t_start = [];
cfg.t_end   = [];

%% ---- T/S quality control ---------------------------------------------
% F2's S stays in [34.33 35.07] over the whole record and T never approaches
% either gate, so the base's [0 33] / [33 36] windows remove NOTHING here.
% They are left on as the guard rather than switched off as the trust -- at
% zero cost to this mooring they would catch a regridding that went wrong.
%
% CHANGE OF BEHAVIOUR WORTH KNOWING: the old f2_config set both to [] and the
% old F2 driver never called qc_ts at all, so this is the first F2 run that
% passes through it. The claim that it removes nothing is the old config's own
% note about F2's salinity range, not something re-measured -- check that the
% first run prints "removed 0 T and 0 S samples". Non-zero means either the
% note was wrong or the product changed, and both are worth knowing.

%% ---- mixed layer -----------------------------------------------------
% The top T sensor is at 78 m and nothing is above it, so the upper 78 m of the
% mean profile is asserted, not measured. 117 m is at the shallowest level that
% carries salinity (see the coverage note below: row 1 at 87 m has no S at all,
% row 2 at 91 m dies in October, row 3 at 107 m is the first usable density),
% which makes the whole genuinely unmeasured cap neutral.
%
% This is the one assumption with real leverage on the answer: c1 spans
% 1.80-2.00 m/s across treatments of the cap, and mode 1's rms p' spikes to
% ~150 Pa at the surface, INSIDE the assumed layer. So this is not a harmless
% boundary condition -- it shapes the quantity the flux is built on.
%
% The FLEAT gliders do not rescue it. See the 22 m depth offset note below:
% anything using glider DEPTH is suspect, so their sigma(z) cannot be spliced
% on. The glider T-S RELATION is fine, being depth-independent, and is the one
% route back into this cap that is still open.
cfg.mld = 50;                  % [m]

%% ---- modes -----------------------------------------------------------
% nbc = 6 from the base. At F2 that is supported: against M2-band eta the
% cumulative variance explained is ~0.86 by four modes and still climbing
% slowly. Against RAW eta it is 0.665 at four modes, which is the WRONG number
% to quote -- raw eta includes the broadband continuum, which four modes were
% never going to represent.

end

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
