function A = align_flux_records(R, dt_hr)
% ALIGN_FLUX_RECORDS  Put several moorings' flux on one common hourly clock.
%
%   A = align_flux_records(R)
%   A = align_flux_records(R, dt_hr)
%
%   R       cell array or struct array of fleat_IT_flux outputs
%   dt_hr   timestep of the common grid [h]. Default 1.
%
%   A .time        [1 x nt]   common grid, datenum, spanning the UNION
%     .name        {1 x nm}
%     .lat .lon .H [1 x nm]
%     .Fx .Fy      [nm x nt]  W/m, NaN where that mooring has no record
%     .mag .dir    [nm x nt]  magnitude [W/m] and direction [deg true, 0 = N]
%     .cover       [nm x 1]   fraction of the common grid each mooring fills
%     .ok          [nm x 1]   logical, which inputs carried a usable flux
%
%   ---------------------------------------------------------------------
%   WHY THE UNION AND NOT THE INTERSECTION
%   ---------------------------------------------------------------------
%   The records do not line up. F2 runs 5 Jun 2016 to 13 Apr 2017; F4 is cut at
%   16 Nov 2016 because its deep ADCP fails. Intersecting would throw away half
%   of F2 to match F4, and would do it silently.
%
%   So the grid spans everything and each mooring is padded with NaN outside its
%   own record. An array-frame movie then shows F4's arrow simply disappearing
%   in November, which is the truth about the data and is visible to a viewer.
%   Anything that averages across moorings has to count the NaNs -- A.cover is
%   there so that a frame with two moorings reporting is not mistaken for a
%   frame with five.
%
%   ---------------------------------------------------------------------
%   WHY ROUNDING AND NOT interp1
%   ---------------------------------------------------------------------
%   These are hourly averages of 10-minute samples, so each record is already
%   on an hour grid, but the records start at different hours and floating-point
%   datenums do not compare exactly. Snapping to the nearest timestep is exact
%   bookkeeping and moves nothing. Interpolating would invent flux values across
%   a gap and across the ends of each record, which is precisely what the NaN
%   padding exists to prevent.
%
%   A mooring whose timestep is NOT dt_hr gets a warning and is snapped anyway:
%   two samples landing in one slot means the later one wins. Nothing in FLEAT
%   is sampled other than hourly, so this is a guard, not a feature.
%
% author TEJ, with help from Claude
% 15 October 2026

if nargin < 2 || isempty(dt_hr), dt_hr = 1; end
if isstruct(R), R = num2cell(R); end
if ~iscell(R) || isempty(R)
    error('align_flux_records:input', ...
        'Pass the cell array of fleat_IT_flux outputs, e.g. align_flux_records(R).');
end

% ---- which inputs are usable ----------------------------------------
nm  = numel(R);
use = false(nm,1);
for i = 1:nm
    r = R{i};
    use(i) = isstruct(r) && isfield(r,'ok') && r.ok && ...
             isfield(r,'flux') && isstruct(r.flux) && ...
             isfield(r,'time') && ~isempty(r.time);
end
if ~any(use)
    error('align_flux_records:none', 'None of the %d inputs carried a flux.', nm);
end

A.name = cell(1,nm);
A.lat  = nan(1,nm);  A.lon = nan(1,nm);  A.H = nan(1,nm);
A.ok   = use;
for i = 1:nm
    if isstruct(R{i}) && isfield(R{i},'name'), A.name{i} = R{i}.name; else, A.name{i} = sprintf('#%d',i); end
    if use(i)
        A.lat(i) = R{i}.lat;  A.lon(i) = R{i}.lon;  A.H(i) = R{i}.H;
    end
end

% ---- the common grid -------------------------------------------------
t1 = inf;  t2 = -inf;
for i = find(use)'
    t = R{i}.time(:).';
    t1 = min(t1, t(1));  t2 = max(t2, t(end));
end
step   = dt_hr/24;
A.time = t1 : step : t2 + step/2;
nt     = numel(A.time);

A.Fx = nan(nm, nt);  A.Fy = nan(nm, nt);
A.cover = zeros(nm,1);

for i = find(use)'
    t  = R{i}.time(:).';
    Fx = R{i}.flux.Fx(:).';
    Fy = R{i}.flux.Fy(:).';

    % the record's own timestep, as a check on the snap below
    if numel(t) > 1
        dt_i = median(diff(t))*24;
        if abs(dt_i - dt_hr) > 1e-6
            warning('align_flux_records:dt', ...
                ['%s is sampled every %.3f h but the common grid is %.3f h. ' ...
                 'Samples are being snapped to the nearest slot and collisions ' ...
                 'keep the later one.'], A.name{i}, dt_i, dt_hr);
        end
    end

    k  = round((t - t1)/step) + 1;
    in = k >= 1 & k <= nt;
    A.Fx(i, k(in)) = Fx(in);
    A.Fy(i, k(in)) = Fy(in);
    A.cover(i) = nnz(isfinite(A.Fx(i,:)))/nt;
end

A.mag = hypot(A.Fx, A.Fy);
A.dir = mod(90 - atan2d(A.Fy, A.Fx), 360);

fprintf('aligned %d of %d moorings on %s to %s (%d steps of %.0f h)\n', ...
        nnz(use), nm, datestr(t1,'dd-mmm-yyyy'), datestr(t2,'dd-mmm-yyyy'), nt, dt_hr);
for i = 1:nm
    if use(i)
        fprintf('  %-3s  %5.1f%% of the grid,  mean |F| %6.3f kW/m\n', ...
                A.name{i}, 100*A.cover(i), mean(A.mag(i,:),'omitnan')/1e3);
    else
        fprintf('  %-3s  no flux -- excluded\n', A.name{i});
    end
end
nboth = sum(isfinite(A.Fx), 1);
fprintf('  moorings reporting simultaneously: %d to %d (median %d)\n', ...
        min(nboth), max(nboth), round(median(nboth)));
end
