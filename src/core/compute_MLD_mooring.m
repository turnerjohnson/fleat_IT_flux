function F_h_out = compute_MLD_mooring(F_h, cfg)
% COMPUTE_MLD_MOORING  Mixed layer depth two ways: density threshold and PEA.
%
%   F_h = compute_MLD_mooring(F_h, cfg)
%
%   F_h   hourly mooring struct AFTER compute_density_mooring, so it carries
%         .depth [nsens x nt] positive down, .sigma (sigma_0) and .time.
%   cfg   mooring config. Every field it reads has a default (see PICK at the
%         bottom), so a config written before this function existed still runs.
%
%   Adds F_h.ml, hourly and [1 x nt] throughout:
%        .z_thresh    density-threshold MLD [m]
%        .z_pea       potential-energy-anomaly MLD [m]
%        .phi         PEA of the layer the PEA method diagnosed [J m^-2], so
%                     cfg.mld_pea_crit by construction where it is resolved
%        .phi_thresh  PEA of the layer the THRESHOLD method diagnosed [J m^-2]
%                     -- the energy that answer implies, Reichl's Figures 3-4
%        .phi_span    PEA over cfg.mld_phi_span metres below the reference
%                     [J m^-2] -- defined whether or not the base is resolved
%        .zref        reference depth actually used that hour [m]
%        .sigma_ref   sigma_0 there [kg m^-3]
%        .dsig_top    sigma_0 jump across the TOP sensor interval [kg m^-3]
%        .unresolved  true when the threshold crossing is inside the top
%                     sensor interval: a bound, not a measurement
%        .pea_pinned  the same for the PEA method -- phi_crit is already spent
%                     one grid step below the reference
%        .nobs        density samples in the profile
%        .info        the settings as run, the counts, and the summary numbers
%
%   ---------------------------------------------------------------------
%   READ THIS FIRST: WHAT A MOORING CAN AND CANNOT SEE
%   ---------------------------------------------------------------------
%   NOTHING MEASURES DENSITY ABOVE THE SHALLOWEST SIGMA SENSOR. At F2 that is
%   107 m (row 1 at 87 m reports no salinity at all); at F4 it is ~100 m. The
%   climatological mixed layer at 8.5 N is 30-80 m, which is ABOVE the array.
%
%   So the reference level here is not 10 m, as in de Boyer Montegut (2004) --
%   it is the top of the measured column, and when the real mixed layer sits
%   above it BOTH methods below return a depth pinned just under the reference
%   and neither has any way to know it. That case is flagged (.unresolved for
%   the threshold, .pea_pinned for the PEA) rather than silently averaged in,
%   and it is the common case here: the F2 field notes
%   record that when the shallowest isopycnal target is missing the surface is
%   above the array 99% of those times.
%
%   What this function is therefore FOR: the hours when the mixed layer does
%   reach below the top sensor, which are the hours that can falsify cfg.mld.
%   cfg.mld is asserted, not measured (see the mixed-layer note in
%   fleat_config_base), it sets the neutral cap in f_background_v2, and c1
%   spans 1.80-2.00 m/s across treatments of that cap. A resolved MLD
%   distribution is the only observational check on it this data can offer.
%   A low resolved fraction is itself the result: it says the cap stays an
%   assertion, and the summary says so in those words.
%
%   An .unresolved hour is an UPPER BOUND, not a measurement: the density is
%   already increasing at the first interval, so the true base is at or above
%   zref. The returned number is kept (it is the bound) and excluded from every
%   statistic printed below.
%
%   ---------------------------------------------------------------------
%   METHOD 1: DENSITY THRESHOLD
%   ---------------------------------------------------------------------
%   MLD is where sigma_0 first reaches sigma_0(zref) + cfg.mld_dsigma, found on
%   the FIRST crossing going down and linearly interpolated within the bracket.
%   dsigma = 0.03 kg/m^3 is de Boyer Montegut's value and is the one number the
%   method has.
%
%   First crossing, not interp1 on the whole column: an hourly profile from 43
%   sensors is not monotonic (Thorpe sorting fires on 29% of F2 profiles), and
%   interp1 against a non-monotonic sigma either errors or picks an arbitrary
%   branch. Within the bracket sigma is increasing by construction, so the
%   interpolation there is well posed.
%
%   ---------------------------------------------------------------------
%   METHOD 2: POTENTIAL ENERGY ANOMALY (Reichl et al. 2022, JGR Oceans,
%   10.1029/2021JC018140 -- in resources/. phi after Simpson 1981.)
%   ---------------------------------------------------------------------
%   phi is the work PER UNIT AREA needed to mix a layer to uniform potential
%   density: its energetic distance from the well-mixed state. Over a layer
%   from zref to h of thickness L = h - zref, with d positive down, Reichl's
%   equation 17 is
%
%       phi(h) = g * int_zref^h (sigma(d) - sigmabar) * d  dd ,
%       sigmabar = (1/L) * int_zref^h sigma dd                   [J m^-2]
%
%   Reichl write it with z measured upward from the surface; d = -z turns it
%   into the form above, which is positive for a stratified layer, zero for a
%   mixed one, and equals g*drho*L^2/12 for a profile linear in sigma_0.
%
%   Four things worth knowing about that expression:
%
%   - IT IS PER AREA, NOT PER VOLUME. Simpson's original phi divides by the
%     layer thickness and is quoted in J/m^3; Reichl's does not, because the
%     point is to compare it with an energy supply -- 10 J/m^2 is a few hours
%     of a 12 m/s wind (their Appendix B). Dividing by L would make the
%     criterion below mean something different at every depth. Keep the units.
%
%   - THE ORIGIN OF THE MOMENT ARM DOES NOT MATTER. int (sigma - sigmabar) dd
%     is zero by construction, so adding any constant to d changes nothing.
%     That is why a layer starting at 107 m rather than at the surface is still
%     a well defined phi, and why the 1000 in rho = 1000 + sigma drops out.
%     Reichl makes the same point about their reference level z_r.
%
%   - IT IS A CUMULATIVE QUANTITY, so it gives a profile phi(h), not one
%     number, and a depth comes from where that profile crosses a criterion.
%
%   - sigma_0 IS THE RIGHT DENSITY HERE, and it is also the approximation.
%     Reichl's exact quantity is dPg, the PE change after homogenising Theta
%     and S separately through the nonlinear equation of state (their eq. 16);
%     phi assumes potential density itself is homogenised (eq. 17). They show
%     the two differ by less than 1% for energies below 1e4 J/m^2 and by more
%     than a few per cent only above 1e5, so phi is the one they recommend:
%     one profile, not two, at no cost in this range. Everything here is well
%     inside it -- the criterion is 10 J/m^2 and .phi_span over 100 m comes out
%     around 1e4 -- and the thermobaric error in sigma_0 that f_background_v2's
%     header documents (it under-reads N^2 by ~24% below 1200 m) is a
%     deep-column effect, small over the layers searched here.
%
%   THE CRITERION. Reichl's proposal is to PRESCRIBE phi and solve for h,
%   rather than prescribe h and report phi:
%
%       MLD_pea = shallowest h where phi(h) >= cfg.mld_pea_crit
%
%   with cfg.mld_pea_crit in J/m^2. 10 is the default here; their sensitivity
%   study recommends 10-25 J/m^2 as the range that reproduces threshold and
%   HT09 MLDs globally, and ~15 for model output. Solving it needs no iteration
%   here because phi(h) is built cumulatively and the crossing is bracketed.
%
%   WHAT THE TWO METHODS DO DIFFERENTLY, exactly. For a profile linear in
%   sigma_0, the two criteria are
%
%       threshold:  h = dsigma / a                 a = dsigma_0/dd [kg m^-4]
%       PEA:        h = (12*phi_crit/(g*a))^(1/3)  since phi = g*a*h^3/12
%
%   which coincide at a = sqrt(g*dsigma^3/(12*phi_crit)), i.e. at h = 20 m for
%   0.03 kg/m^3 and 10 J/m^2. DEEPER than that the PEA criterion is the
%   stricter of the two and returns a shallower MLD; shallower than that it is
%   looser. That asymmetry is Reichl's whole argument: because the energy cost
%   of mixing grows with the DISTANCE fluid is moved as well as with the
%   density contrast, a fixed dsigma admits layers whose phi ranges over
%   70-370 J/m^2 (their Figure 1) and exceeds 1500 J/m^2 for deep winter
%   layers, so the threshold method calls columns "mixed" that no plausible
%   energy supply could have mixed. The threshold MLD is a statement about two
%   neighbouring samples; the PEA MLD is a statement about the whole layer.
%
%   Where the two agree, the base is a sharp step and the number is robust.
%   Where they disagree, the profile has no step, "the" MLD is genuinely
%   ambiguous, and .phi_thresh below says how much energy the threshold answer
%   was implicitly asserting. Report both.
%
%   The alternative criterion, cfg.mld_pea_crit = [], drops Reichl's fixed
%   energy for the LINEAR EQUIVALENT of the same dsigma the threshold uses,
%   g*dsigma*L^2/12, evaluated at each h. That makes both methods run off one
%   physical parameter, so their difference isolates the SHAPE of the profile
%   with no second constant involved -- useful for the comparison, not for
%   quoting an MLD, since it inherits exactly the energetic inconsistency
%   Reichl is objecting to.
%
%   ---------------------------------------------------------------------
%   WHY phi IS BUILT ON A FINE GRID
%   ---------------------------------------------------------------------
%   phi(h) is evaluated on a cfg.mld_dz grid between the samples, not only at
%   the sensors. The interpolation is LINEAR, which is exactly the profile the
%   trapezoidal integral already assumes, so this adds no information and no
%   smoothing -- it only resolves where the crossing falls. At the sensors
%   alone, with 10-148 m spacing, the crossing depth would be quantised to the
%   sensor grid and the comparison with the interpolated threshold depth would
%   be meaningless.
%
% as part of pipeline for computing IT flux on FLEAT moorings
%
% author TEJ, with help from Claude
% 22 October 2026

%% ---- config, with defaults -------------------------------------------
g      = 9.81;                                  % [m s^-2]
dsig   = pick(cfg, 'mld_dsigma',  0.03);        % [kg m^-3]

% NOT pick(), on purpose. For this one field [] is a MEANING -- use the linear
% equivalent of dsigma instead of a fixed energy -- and pick treats an empty
% field as absent and would hand back the default instead.
if isstruct(cfg) && isfield(cfg, 'mld_pea_crit')
    phicr = cfg.mld_pea_crit;                   % [J m^-2], or [] (see header)
else
    phicr = 10;                                 % Reichl et al. 2022
end
zrefin = pick(cfg, 'mld_zref',     []);         % [] = shallowest sample, hourly
zmax   = pick(cfg, 'mld_zmax',    400);         % [m] below zref, search limit
span   = pick(cfg, 'mld_phi_span', 100);        % [m] layer for the phi index
dzg    = pick(cfg, 'mld_dz',        1);         % [m] grid for phi(h)
minobs = pick(cfg, 'mld_minobs',    4);         % density samples in a profile
mname  = pick(cfg, 'mooring', 'mooring');

if ~(isscalar(dsig) && isfinite(dsig) && dsig > 0)
    error('compute_MLD_mooring:dsigma', ...
        'cfg.mld_dsigma must be a positive scalar [kg m^-3]; got %s.', mat2str(dsig));
end

for q = {'depth','sigma','time'}
    if ~isfield(F_h, q{1}) || isempty(F_h.(q{1}))
        error('compute_MLD_mooring:input', ...
            ['F_h.%s is missing. This runs on the hourly struct AFTER ' ...
             'compute_density_mooring, which is what puts .sigma there.'], q{1});
    end
end
D = F_h.depth;  SG = F_h.sigma;  t = F_h.time(:).';
if ~isequal(size(D), size(SG))
    error('compute_MLD_mooring:size', ...
        'F_h.depth is %s and F_h.sigma is %s; they must match.', ...
        mat2str(size(D)), mat2str(size(SG)));
end
nt = numel(t);
if size(D,2) ~= nt
    error('compute_MLD_mooring:clock', ...
        'F_h.depth has %d columns against %d times.', size(D,2), nt);
end
if ~isreal(SG)
    error('compute_MLD_mooring:complex', ...
        ['F_h.sigma is complex, so every comparison below is meaningless. ' ...
         'See the salinity note in compute_density_mooring.']);
end

%% ---- hour by hour ----------------------------------------------------
ML = struct('time', t, ...
            'z_thresh',   nan(1,nt), 'z_pea', nan(1,nt), 'phi', nan(1,nt), ...
            'phi_thresh', nan(1,nt), 'phi_span', nan(1,nt), ...
            'zref',       nan(1,nt), 'sigma_ref', nan(1,nt), ...
            'dsig_top',   nan(1,nt), 'unresolved', false(1,nt), ...
            'pea_pinned', false(1,nt), 'nobs', zeros(1,nt));
n = struct('short', 0, 'noref', 0, 'thin', 0, 'nocross_t', 0, 'nocross_p', 0);

for it = 1:nt
    ok = isfinite(D(:,it)) & isfinite(SG(:,it));
    if nnz(ok) < minobs
        n.short = n.short + 1;  continue
    end

    % Sorted, one sample per depth. Duplicate depths are not expected -- the
    % product interpolates each sensor's own depth -- and if two coincide the
    % shallower-listed one wins rather than breaking the interpolation.
    sv      = SG(ok,it);
    [d, iu] = unique(D(ok,it));
    s       = sv(iu);
    if numel(d) < 2
        n.thin = n.thin + 1;  continue     % one depth: nothing to integrate
    end

    % THE REFERENCE. Default is the shallowest sample of THAT hour. Knockdown
    % moves it (F2's sensors range 56 m in the median, 63 m at most) and so does
    % sensor dropout, so .zref is returned per hour and its spread reported.
    % A fixed
    % cfg.mld_zref is the alternative, and it costs every hour whose top sample
    % has been pushed below that depth -- which at F2 is most of a knockdown
    % event, so it is not the default.
    if isempty(zrefin)
        zr = d(1);  sr = s(1);
    else
        if d(1) > zrefin || d(end) < zrefin
            n.noref = n.noref + 1;  continue
        end
        zr   = zrefin;
        sr   = interp1(d, s, zr);
        keep = d > zr;
        d    = [zr; d(keep)];
        s    = [sr; s(keep)];
    end
    if numel(d) < 2 || ~isfinite(sr)
        n.thin = n.thin + 1;  continue
    end

    ML.zref(it)      = zr;
    ML.sigma_ref(it) = sr;
    ML.nobs(it)      = numel(d);
    ML.dsig_top(it)  = s(2) - s(1);

    % ---- method 1: density threshold ---------------------------------
    tgt = sr + dsig;
    k   = find(s >= tgt, 1);            % k >= 2 always, since s(1) = sr < tgt
    if isempty(k)
        n.nocross_t = n.nocross_t + 1;
    else
        ML.z_thresh(it) = d(k-1) + (tgt - s(k-1)) * (d(k) - d(k-1)) / (s(k) - s(k-1));
        % The crossing is inside the top sensor interval: nothing above it
        % brackets a uniform layer, so the base is at or above zref and this
        % depth is an upper bound. Also true whenever dsig_top >= dsigma.
        ML.unresolved(it) = (k == 2);
    end

    % ---- method 2: potential energy anomaly --------------------------
    zb = min(d(end), zr + zmax);
    zg = (zr : dzg : zb)';
    if numel(zg) < 3
        n.thin = n.thin + 1;  continue     % less than 2 m of column below zref
    end
    sg = interp1(d, s, zg);                     % the profile trapz assumes anyway

    L   = zg - zr;
    I1  = cumtrapz(zg, sg);                     % int sigma dd
    I2  = cumtrapz(zg, sg .* zg);               % int sigma*d dd
    Dm  = (zg.^2 - zr^2) / 2;                   % int d dd, exactly
    phi = nan(size(zg));
    j   = L > 0;
    % Reichl eq. 17, per unit AREA. sigmabar = I1/L is the layer mean, and the
    % moment arm needs no origin shift because int (sigma - sigmabar) dd = 0.
    phi(j) = g * (I2(j) - (I1(j)./L(j)) .* Dm(j));

    % The stratification index that survives the censoring. phi over a FIXED
    % thickness below the reference needs no mixed layer base to be resolved,
    % so it is defined in every hour that has a profile, and it is comparable
    % across hours and across moorings in a way the MLD here is not.
    if zg(end) >= zr + span
        ML.phi_span(it) = interp1(zg(j), phi(j), zr + span);
    end

    % HOW MUCH ENERGY THE THRESHOLD ANSWER WAS ASSERTING. This is the quantity
    % Reichl's Figures 3 and 4 map: phi evaluated at the threshold MLD. Its
    % spread across hours is the criticism, in one number.
    if isfinite(ML.z_thresh(it)) && ML.z_thresh(it) <= zg(end)
        ML.phi_thresh(it) = interp1(zg(j), phi(j), max(ML.z_thresh(it), zg(2)));
    end

    if isempty(phicr)
        crit = g * dsig * L.^2 / 12;            % linear-profile equivalent
    else
        crit = phicr * ones(size(L));
    end

    f    = phi - crit;
    f(1) = -Inf;                                % L = 0: no layer, no energy
    kk   = find(f >= 0, 1);
    if isempty(kk)
        n.nocross_p = n.nocross_p + 1;
    elseif kk == 2
        % One grid step below the reference already costs phi_crit, so this is
        % the PEA method's version of .unresolved: a bound, not a measurement.
        ML.z_pea(it)      = zg(2);
        ML.phi(it)        = phi(2);
        ML.pea_pinned(it) = true;
    else
        w = f(kk-1) / (f(kk-1) - f(kk));
        ML.z_pea(it) = zg(kk-1) + w*(zg(kk) - zg(kk-1));
        ML.phi(it)   = phi(kk-1) + w*(phi(kk) - phi(kk-1));
    end
end

%% ---- what the record says --------------------------------------------
% Every statistic here is over hours RESOLVED by the method in question. An
% unresolved hour is a bound, and averaging bounds in with measurements would
% manufacture a mixed layer depth out of the depth of the top sensor.
res_t = isfinite(ML.z_thresh) & ~ML.unresolved;
res_p = isfinite(ML.z_pea)    & ~ML.unresolved & ~ML.pea_pinned;
both  = res_t & res_p;

I = struct();
I.dsigma      = dsig;
I.pea_crit    = phicr;              % [] when the linear equivalent was used
I.zref_fixed  = zrefin;
I.zmax        = zmax;
I.phi_span_m  = span;
I.phi_span_med = median(ML.phi_span, 'omitnan');
I.dz          = dzg;
I.minobs      = minobs;
I.counts      = n;
I.nt          = nt;
I.frac_res_t  = mean(res_t);
I.frac_res_p  = mean(res_p);
I.zref_med    = median(ML.zref, 'omitnan');
I.zref_range  = [min(ML.zref) max(ML.zref)];
I.z_thresh_med = median(ML.z_thresh(res_t), 'omitnan');
I.z_pea_med    = median(ML.z_pea(res_p),    'omitnan');
I.z_thresh_pct = prct(ML.z_thresh(res_t), [10 90]);
I.z_pea_pct    = prct(ML.z_pea(res_p),    [10 90]);
% Reichl's criticism, as a number for this record: the energy the threshold
% MLD implies, and how far it ranges. A tight distribution would mean the two
% methods are interchangeable here; a wide one means they are not.
I.phi_thresh_med = median(ML.phi_thresh(res_t), 'omitnan');
I.phi_thresh_pct = prct(ML.phi_thresh(res_t), [10 90]);
I.pea_minus_thresh = median(ML.z_pea(both) - ML.z_thresh(both), 'omitnan');
I.ratio_med        = median(ML.z_pea(both) ./ ML.z_thresh(both), 'omitnan');
if nnz(both) > 2
    cc = corrcoef(ML.z_thresh(both), ML.z_pea(both));
    I.corr = cc(1,2);
else
    I.corr = NaN;
end
ML.info = I;

fprintf('[%s] MLD: reference %.0f m median (%.0f to %.0f m over the record)\n', ...
        mname, I.zref_med, I.zref_range(1), I.zref_range(2));
fprintf('       resolved: threshold %.0f%% of %d hours, PEA %.0f%%, both %.0f%%\n', ...
        100*I.frac_res_t, nt, 100*I.frac_res_p, 100*mean(both));
if any(both)
    fprintf('       threshold (dsigma %.3f): median %.0f m, 10-90%% %.0f-%.0f m\n', ...
            dsig, I.z_thresh_med, I.z_thresh_pct(1), I.z_thresh_pct(2));
    if isempty(phicr)
        fprintf('       PEA (phi_crit = g*%.3f*L^2/12): median %.0f m, 10-90%% %.0f-%.0f m\n', ...
                dsig, I.z_pea_med, I.z_pea_pct(1), I.z_pea_pct(2));
    else
        fprintf('       PEA (phi_crit = %.3g J m^-2, Reichl 2022): median %.0f m, 10-90%% %.0f-%.0f m\n', ...
                phicr, I.z_pea_med, I.z_pea_pct(1), I.z_pea_pct(2));
    end
    fprintf('       PEA minus threshold: median %+.0f m (ratio %.2f, r = %.2f)\n', ...
            I.pea_minus_thresh, I.ratio_med, I.corr);
    fprintf(['       phi AT the threshold MLD: median %.0f J m^-2, 10-90%% %.0f-%.0f ' ...
             '(Reichl Fig 3-4: a wide range means the two are not interchangeable)\n'], ...
             I.phi_thresh_med, I.phi_thresh_pct(1), I.phi_thresh_pct(2));
end
fprintf('       phi over %.0f m below the reference: median %.3g J m^-2 (%d%% of hours)\n', ...
        span, I.phi_span_med, round(100*mean(isfinite(ML.phi_span))));
if n.short + n.noref + n.thin > 0
    fprintf('       skipped: %d hour(s) with < %d density samples', n.short, minobs);
    if n.noref > 0, fprintf(', %d with no sample at the fixed reference', n.noref); end
    if n.thin  > 0, fprintf(', %d with nothing below the reference', n.thin); end
    fprintf('\n');
end
if n.nocross_t > 0 || n.nocross_p > 0
    fprintf(['       no crossing within %.0f m of the reference: %d hour(s) ' ...
             'threshold, %d PEA (a column mixed deeper than the search)\n'], ...
             zmax, n.nocross_t, n.nocross_p);
end

%% ---- against the cap the modes actually use --------------------------
% This is the comparison the function exists for. cfg.mld is the depth above
% which f_background_v2 sets N^2 to cfg.N2_ml, and it is asserted from the
% shallowest sigma-bearing sensor rather than measured.
if isfield(cfg,'mld') && ~isempty(cfg.mld) && isfinite(cfg.mld) && any(res_t)
    deeper = mean(ML.z_thresh(res_t) > cfg.mld);
    fprintf(['       cfg.mld = %.0f m; the measured base is deeper than that in ' ...
             '%.0f%% of resolved hours\n'], cfg.mld, 100*deeper);
end
if I.frac_res_t < 0.10
    warning('compute_MLD_mooring:unresolved', ...
        ['%s resolves a mixed layer base in only %.1f%% of hours: the ' ...
         'shallowest density is at %.0f m and the layer is above it. cfg.mld ' ...
         '(%.0f m) therefore remains an ASSERTION -- this run does not check ' ...
         'it, and the c1 spread across treatments of the cap stands as the ' ...
         'uncertainty. Do not quote a mean MLD from this record.'], ...
         mname, 100*I.frac_res_t, I.zref_med, pick(cfg,'mld',NaN));
end

F_h.ml  = ML;
F_h_out = F_h;
end

%% ======================================================================
function v = pick(cfg, name, default)
% PICK  cfg field if it exists and is non-empty, otherwise the default. Same
%       helper as f_background_v2, for the same reason: this function must run
%       against the five configs as they stand today.
if isstruct(cfg) && isfield(cfg, name) && ~isempty(cfg.(name))
    v = cfg.(name);
else
    v = default;
end
end

%% ======================================================================
function p = prct(x, q)
% PRCTILE without the Statistics Toolbox, on the finite values only. Returns
% NaN for an empty sample rather than erroring, so the summary still prints.
x = sort(x(isfinite(x)));
if isempty(x), p = nan(size(q)); return; end
if isscalar(x), p = x*ones(size(q)); return; end
% Linear interpolation on the (i-0.5)/n plotting positions, which is MATLAB's
% own prctile convention, so the numbers match if anyone checks them there.
pp = 100*((1:numel(x)) - 0.5)/numel(x);
p  = interp1(pp, x, q, 'linear');
p(q < pp(1))   = x(1);
p(q > pp(end)) = x(end);
end
