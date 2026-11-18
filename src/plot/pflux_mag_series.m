function h = pflux_mag_series(R, opt)
%PFLUX_MAG_SERIES  |F| for every mooring on one common clock: the movie's
%                  lower panel, static, with a running mean over it.
%
%   h = pflux_mag_series()                  % all five from the cache
%   h = pflux_mag_series(R)                 % R as passed to array_map_movie
%   h = pflux_mag_series({'F2','F3'}, opt)
%
%   R    cell/struct array of fleat_IT_flux outputs, OR a cellstr of mooring
%        names to load from data/flux_<name>.mat, OR [] for all five.
%   opt .smooth_days  the heavy line, a centred running mean [d]. Default 3.
%       .min_cover    fraction of a smoothing window that must be finite for
%                     the heavy line to be drawn there. Default 0.6. See
%                     WHAT THE HEAVY LINE IS NOT ALLOWED TO DO.
%       .t0 .t1       window, datenum. [] [] is the whole union record.
%       .ruler        14.765-d verticals. Default true. Phase is ARBITRARY.
%       .kW           kW/m rather than W/m. Default true.
%       .logy         log y-axis. Default false. Try it when one mooring is an
%                     order of magnitude above the rest and the fortnight in
%                     the weak ones is a flat line at the bottom.
%       .dt_hr        common clock. Default 1, i.e. what the records are.
%       .cachedir     where flux_<name>.mat lives. Default ../../data.
%       .fignum       default 22.
%
%   h .fig .ax .A .mag .sm .names .peaks
%
%   ---------------------------------------------------------------------
%   WHAT THIS IS
%   ---------------------------------------------------------------------
%   Exactly the series panel of ARRAY_MAP_MOVIE -- same ALIGN_FLUX_RECORDS
%   clock, same |F| = hypot(Fx, Fy), same one-colour-per-mooring -- drawn as a
%   figure in its own right, at full height, with the record end to end and
%   the running mean the movie has no room for. The movie's panel exists to
%   carry a cursor; this one exists to be read.
%
%   thin   |F| off the pipeline: the 25-h (2 M2 period) phase average that
%          COMPUTE_TIDAL_FLUX already did. It is not a smoothing this figure
%          chose, and it is the smallest average for which the word "flux"
%          means anything (see WHY THE AVERAGE IS NOT OPTIONAL there).
%   heavy  a further opt.smooth_days running mean of that, which is the line
%          the fortnight should be visible in.
%
%   ---------------------------------------------------------------------
%   WHY THE MEAN IS OF |F| AND NOT OF THE VECTOR
%   ---------------------------------------------------------------------
%   mean|F| and |mean F| answer different questions. This figure smooths the
%   MAGNITUDE, so the heavy line is the typical strength of the flux over
%   three days regardless of where it pointed. Smoothing Fx and Fy first and
%   taking the magnitude after would let three days of rotating flux cancel to
%   near zero, and a dip in that line would mean "the direction turned", not
%   "the wave got weaker". Since the question here is whether the STRENGTH is
%   modulated fortnightly, the magnitude is the thing to average.
%
%   The consequence is that this line cannot go to zero and its mean is biased
%   high: |F| is a positive quantity, so noise in Fx and Fy raises it. Read
%   peaks against troughs, not against zero.
%
%   ---------------------------------------------------------------------
%   WHAT THE HEAVY LINE IS NOT ALLOWED TO DO
%   ---------------------------------------------------------------------
%   movmean(..., 'omitnan') paints straight across a gap: a 3-day window with
%   four good hours in it returns the mean of those four hours and looks like
%   every other point on the line. F4 stops in November and every record has
%   hours where the fit failed, so the line would run on through holes it has
%   no business in.
%
%   So the window's finite fraction is counted and the heavy line is cut
%   wherever it falls below opt.min_cover. The thin line is left exactly as
%   the pipeline produced it, gaps and all, so the two together show both the
%   smooth trend and where it is and is not supported.
%
%   ---------------------------------------------------------------------
%   THE 14.765-DAY RULER, AND WHAT COUNTS AS AN ANSWER
%   ---------------------------------------------------------------------
%   Spring/neap is the M2-S2 beat, 1/(2 - 1.9323) cpd = 14.765 d. The dotted
%   verticals are spaced by it and ANCHORED AT THE WINDOW START, not at
%   syzygy: nothing in this repository knows the lunar phase. They measure
%   SPACING only. Peaks that keep step with the ticks across the record are
%   fortnightly; peaks that drift steadily across them are not, and three
%   peaks at roughly a fortnight is the least that distinguishes a cycle from
%   a run of events.
%
%   The console lists the local maxima of each heavy line and the intervals
%   between them, which is the same reading done with a number instead of an
%   eye. A mooring whose intervals sit near 14.8 d has a spring/neap signal in
%   |F|; one whose intervals scatter from 6 to 25 d does not, however
%   rhythmic the line looks.
%
% author TEJ, with help from Claude
% 12 November 2026

%% ---- arguments -------------------------------------------------------
if nargin < 1, R   = []; end
if nargin < 2, opt = struct(); end

def = struct('smooth_days',3, 'min_cover',0.6, 't0',[], 't1',[], ...
             'ruler',true, 'kW',true, 'logy',false, 'dt_hr',1, ...
             'cachedir','', 'fignum',22);
for f = fieldnames(def)'
    % isempty(false) is false, so a logical option set to false survives this.
    if ~isfield(opt, f{1}) || isempty(opt.(f{1})), opt.(f{1}) = def.(f{1}); end
end
if isempty(opt.cachedir)
    opt.cachedir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data');
end

if isempty(R), R = {'F1','F2','F3','F4','F6'}; end
if isstruct(R), R = num2cell(R); end

% A cellstr of names means load; a cell of structs is already loaded.
if iscellstr(R) %#ok<ISCLSTR>
    names = R;  R = cell(1, numel(names));
    for i = 1:numel(names)
        fn = fullfile(opt.cachedir, sprintf('flux_%s.mat', names{i}));
        if ~isfile(fn)
            error('pflux_mag_series:nocache', ...
                ['No %s. Run turner_fleat_IT_flux section 1 for %s first, ' ...
                 'or pass the in-memory R.'], fn, names{i});
        end
        R{i} = load(fn).out;
        fprintf('[%s] from %s\n', names{i}, fn);
    end
end

%% ---- the common clock ------------------------------------------------
% ALIGN_FLUX_RECORDS owns the union grid, the snap-not-interpolate rule and
% the NaN padding. Going through it is what makes this figure the movie's
% panel rather than a second alignment that could disagree with it.
A  = align_flux_records(R, opt.dt_hr);
nm = numel(A.name);

sc = 1;  ulab = 'W m^{-1}';
if opt.kW, sc = 1e-3;  ulab = 'kW m^{-1}'; end
mag = A.mag * sc;                      % [nm x nt], NaN outside each record

if isempty(opt.t0), opt.t0 = A.time(1);   end
if isempty(opt.t1), opt.t1 = A.time(end); end
if opt.t1 <= opt.t0
    error('pflux_mag_series:window', ...
        'opt.t1 (%s) is not after opt.t0 (%s).', ...
        datestr(opt.t1,'dd-mmm-yyyy'), datestr(opt.t0,'dd-mmm-yyyy'));
end
ncyc = (opt.t1 - opt.t0)/14.765;

%% ---- the running mean -----------------------------------------------
% Centred, in samples of the common clock, with the window's finite fraction
% counted so the line can be cut where it is unsupported. See the header.
nw = max(3, round(opt.smooth_days*24/opt.dt_hr));
sm = nan(size(mag));
if opt.smooth_days > 0
    fin  = isfinite(mag);
    frac = movmean(double(fin), nw, 2, 'Endpoints','shrink');
    sm   = movmean(mag, nw, 2, 'omitnan', 'Endpoints','shrink');
    sm(frac < opt.min_cover | ~any(fin,2)) = NaN;
end
fprintf(['running mean: %g d = %d samples of %g h, cut where under %.0f%% of ' ...
         'the window is finite\n'], opt.smooth_days, nw, opt.dt_hr, 100*opt.min_cover);

%% ---- draw ------------------------------------------------------------
% lines(), in the same order as A.name, so a mooring is the colour it is in
% ARRAY_MAP_MOVIE's panel and in PFLUX_UVP_SERIES.
cM = lines(max(nm,7));

h.fig = figure(opt.fignum); clf
set(h.fig, 'Color','w', 'Position',[80 140 1060 520])
ax = axes(h.fig, 'Position',[0.08 0.13 0.88 0.76]);
hold(ax,'on'); box(ax,'on')

% the ruler first, so the data sits on top of it
if opt.ruler
    for tk = opt.t0 : 14.765 : opt.t1
        xline(ax, tk, ':', 'Color',[.6 .6 .6], 'LineWidth',0.8);
    end
end

% The legend names EVERY line, not one per mooring: the two traces are two
% processings of the same series and which is which is the thing a reader has
% to know first. That belongs beside the line, not in a subtitle.
hl  = gobjects(1,0);
lab = cell(1,0);
for i = 1:nm
    if ~any(isfinite(mag(i,:)))
        hl(end+1)  = plot(ax, NaN, NaN, '-', 'Color',cM(i,:), 'LineWidth',1.8); %#ok<AGROW>
        lab{end+1} = sprintf('%s  (no flux)', A.name{i}); %#ok<AGROW>
        continue
    end

    % Thin: the pipeline's own 25-h phase average, at 40% alpha so five of
    % them do not become a solid band the heavy lines are lost in.
    thin = plot(ax, A.time, mag(i,:), '-', 'Color',[cM(i,:) 0.40], 'LineWidth',0.8);

    % The boxcar compute_tidal_flux actually used, in hours, read back off cfg
    % rather than hard-coded: 2 periods at hourly sampling is 25 h, but
    % cfg.flux_periods is a knob and a label that lies is worse than none.
    nper = 25;
    if isfield(R{i},'cfg') && all(isfield(R{i}.cfg, {'M2_cpd','dt_hr'}))
        per = 2;
        if isfield(R{i}.cfg,'flux_periods') && ~isempty(R{i}.cfg.flux_periods)
            per = R{i}.cfg.flux_periods;
        end
        nper = round(per * (24/R{i}.cfg.M2_cpd) / R{i}.cfg.dt_hr) * R{i}.cfg.dt_hr;
    end
    tlab = sprintf('%s  %g-h phase average', A.name{i}, nper);
    if A.cover(i) < 0.995
        tlab = sprintf('%s  (%.0f%% of the grid)', tlab, 100*A.cover(i));
    end
    hl(end+1) = thin;  lab{end+1} = tlab; %#ok<AGROW>

    if opt.smooth_days > 0 && any(isfinite(sm(i,:)))
        heav = plot(ax, A.time, sm(i,:), '-', 'Color',cM(i,:), 'LineWidth',2.0);
        hl(end+1)  = heav; %#ok<AGROW>
        lab{end+1} = sprintf('%s  %g-d running mean of |F|', ...
                             A.name{i}, opt.smooth_days); %#ok<AGROW>
    end
end

xlim(ax, [opt.t0 opt.t1])
if opt.logy
    set(ax, 'YScale','log');
else
    % |F| >= 0, so a floor below zero is empty paper. Guarded: an all-NaN
    % panel leaves yl(2) at the default and [0 0] is not a valid limit.
    yl = ylim(ax);
    if yl(2) > 0, ylim(ax, [0 yl(2)]); end
end
grid(ax,'on')
ylabel(ax, sprintf('|F| [%s]', ulab))
if opt.t1 - opt.t0 <= 75
    datetick(ax, 'x', 'dd mmm', 'keeplimits');
else
    datetick(ax, 'x', 'mmm yyyy', 'keeplimits');
end
xlabel(ax, sprintf('%s to %s', datestr(opt.t0,'dd mmm yyyy'), ...
                               datestr(opt.t1,'dd mmm yyyy')))
% Docked under the axes, not inside them: two entries per mooring is ten
% entries for the full array and no corner of this panel is empty for that
% long. Outside the axes it costs height and never covers the series.
h.lgd = legend(ax, hl, lab, 'Box','off', 'FontSize',8.5, ...
               'NumColumns', max(nm,1), 'Location','southoutside');
title(ax, {'magnitude of the depth-integrated baroclinic M_2 energy flux', ...
           sprintf('%s  |  %.0f d = %.1f spring/neap cycles%s', ...
                   strjoin(A.name(A.ok(:)'), ', '), opt.t1-opt.t0, ncyc, ...
                   tern(opt.ruler, '  |  dotted: 14.765 d, phase arbitrary', ''))}, ...
      'FontWeight','bold')

h.ax = ax;  h.A = A;  h.mag = mag;  h.sm = sm;  h.names = A.name;
h.units = ulab;

%% ---- the reading, as numbers -----------------------------------------
% The eye on the heavy line, done with peak spacing. A cycle is a SPACING
% claim, so the intervals are what gets printed, not the peak heights.
fprintf('\n%s to %s  (%.0f d, %.1f spring/neap cycles of 14.765 d)\n', ...
        datestr(opt.t0,'dd-mmm-yyyy'), datestr(opt.t1,'dd-mmm-yyyy'), ...
        opt.t1-opt.t0, ncyc);
in = A.time >= opt.t0 - 1e-6 & A.time <= opt.t1 + 1e-6;
h.peaks = cell(1,nm);
for i = 1:nm
    g = mag(i,in);
    if ~any(isfinite(g))
        fprintf('  %-3s  nothing in this window\n', A.name{i});  continue
    end
    gf = g(isfinite(g));
    fprintf('  %-3s  mean %6.3f  rms %6.3f  max %6.3f %s\n', ...
            A.name{i}, mean(gf), sqrt(mean(gf.^2)), max(gf), ulab);

    s  = sm(i,in);  t = A.time(in);
    tp = local_max(t, s, 5);                % 5 d minimum separation
    h.peaks{i} = tp;
    if numel(tp) < 2
        fprintf('       %d peak(s) in the %g-d mean: too few to space\n', ...
                numel(tp), opt.smooth_days);
        continue
    end
    d = diff(tp);
    fprintf('       %d peaks, intervals %s d  (mean %.1f, fortnight is 14.8)\n', ...
            numel(tp), strjoin(compose('%.1f', d(:)'), ' '), mean(d));
    % Modulation depth of the smooth line, which is what "is there a cycle"
    % amounts to if the spacing checks out: a 10% swing is not a cycle.
    sf = s(isfinite(s));
    fprintf('       %g-d mean runs %.3f to %.3f %s (peak/trough %.1f)\n', ...
            opt.smooth_days, min(sf), max(sf), ulab, max(sf)/max(min(sf),eps));
end
fprintf(['\nRead the heavy lines: peaks a fortnight apart, keeping step with the\n' ...
         'dotted ticks, are spring/neap. Peaks that drift across them are not.\n\n']);
end

%% ======================================================================
function tp = local_max(t, y, minsep_d)
%LOCAL_MAX  Times of the local maxima of y(t), no two closer than minsep_d.
%
%   Strongest first, each one blanking its neighbourhood, which is the same
%   rule findpeaks('MinPeakDistance') uses -- written out because findpeaks is
%   in the Signal Processing Toolbox and this figure should not need it. NaN
%   runs break the series, so a peak is never reported across a gap.
tp = [];
if numel(y) < 3, return, end
y   = y(:)';  t = t(:)';
cand = find(isfinite(y(2:end-1)) & isfinite(y(1:end-2)) & isfinite(y(3:end)) ...
            & y(2:end-1) > y(1:end-2) & y(2:end-1) >= y(3:end)) + 1;
if isempty(cand), return, end
[~, ord] = sort(y(cand), 'descend');
cand = cand(ord);
keep = [];
for k = cand
    if all(abs(t(k) - t(keep)) >= minsep_d), keep(end+1) = k; end %#ok<AGROW>
end
tp = sort(t(keep));
end

%% ======================================================================
function y = tern(c, a, b)
if c, y = a; else, y = b; end
end
