function h = pflux_uvp_series(R, t0, opt)
%PFLUX_UVP_SERIES  <u'p'> and <v'p'> through time, one panel each, by mooring.
%
%   h = pflux_uvp_series()                       % F2 and F3 from the cache
%   h = pflux_uvp_series({'F2','F3'}, t0)
%   h = pflux_uvp_series(R, t0, opt)             % R from turner_fleat_IT_flux
%
%   R    cell/struct array of fleat_IT_flux outputs, OR a cellstr of mooring
%        names to load from data/flux_<name>.mat, OR [] for {'F2','F3'}.
%   t0   window start, datenum. [] starts at the first hour every mooring has.
%   opt .ndays        window length [d]. Default 45 -- see WHY 45 DAYS.
%       .t_end        datenum, overrides ndays.
%       .smooth_days  the heavy overlay line, a running mean [d]. Default 3.
%                     0 draws only the 25-h phase-averaged flux.
%       .raw          true (default) scatters the instantaneous product
%                     underneath, where the record carries it. See THE THREE
%                     LINES before reading anything off the dots.
%       .ruler        true (default) marks 14.765-d intervals from t0.
%       .sameY        true (default) gives both panels one y-axis.
%       .kW           true (default) plots kW/m rather than W/m.
%       .cachedir     where flux_<name>.mat lives. Default ../../data.
%       .fignum       default 21.
%
%   h .fig .ax(2) .lgd .names .t0 .t1
%
%   ---------------------------------------------------------------------
%   THE LAYOUT
%   ---------------------------------------------------------------------
%   Two stacked panels, ONE PER COMPONENT: <u'p'> on top, <v'p'> below, with
%   every mooring drawn in both and COLOUR CARRYING THE MOORING. So a panel
%   answers "do these sites agree?" and reading down the stack answers "which
%   way is it going?".
%
%   Both panels share one y-axis by default (opt.sameY). u and v are the same
%   quantity resolved on two axes and giving them separate rulers would make a
%   weak cross-shore flux look like a strong one.
%
%   ---------------------------------------------------------------------
%   WHAT IS ON THE AXES
%   ---------------------------------------------------------------------
%   out.flux.Fx and .Fy, which are int u'p' dz and int v'p' dz -- the
%   DEPTH-INTEGRATED baroclinic M2 flux in W/m, already phase-averaged by
%   compute_tidal_flux over 2 M2 periods (a 25-h boxcar). So the 2*omega
%   wiggle is gone before it reaches this figure and what is left is the
%   slowly varying flux: exactly the quantity whose fortnightly envelope is
%   the thing being looked for. Plotting the instantaneous u'p' instead would
%   bury the fortnight under a 6-hourly oscillation twice its size.
%
%   Sign is geographic, not along-beam: +Fx is eastward, +Fy northward. A
%   series that crosses zero is a flux that reverses, which is a standing
%   wave or a change of source, and it matters more than the magnitude.
%
%   ---------------------------------------------------------------------
%   THE THREE LINES, AND WHY THE DOTS ARE NOT DATA
%   ---------------------------------------------------------------------
%   dots   int u'p' dz with NO phase average: flux.Fx_inst, the thing
%          compute_tidal_flux averages. Faint, behind everything.
%   thin   flux.Fx itself -- the 25-h (2 M2 period) boxcar. This is the
%          pipeline's output and the figure's subject; it is not a smoothing
%          this plot chose.
%   heavy  a further opt.smooth_days running mean, for the eye only.
%
%   The dots are the raw PRODUCT of two fitted amplitudes, not a measurement.
%   For two narrowband M2 signals their product is
%
%       A cos(wt + phi_u) cos(wt + phi_p)
%           = (A/2)[ cos(phi_u - phi_p) + cos(2wt + phi_u + phi_p) ]
%
%   so the dots oscillate about the thin line at TWICE the tidal frequency
%   with an amplitude of order the flux itself, and spend half of every tidal
%   cycle on the wrong side of zero. An instantaneous dot of -2 kW/m at a site
%   whose flux is +1 kW/m east is not a westward flux; it is the cos(2wt) term
%   at its trough. Only the average over whole periods is a flux.
%
%   They are worth drawing for two things the thin line cannot show: whether
%   the scatter is symmetric about the average (it should be, if the signal
%   really is narrowband M2 and the boxcar is doing its job), and where the
%   record is actually sampled, since a hour with no fit leaves a gap in the
%   dots that the boxcar paints straight over.
%
%   WHERE THEY COME FROM. flux.Fx_inst was added to compute_tidal_flux on
%   12 Nov 2026, so a flux_<name>.mat cached before that has no dots: rerun
%   the mooring, or pass an in-memory out from opt.keep = 'full' and this
%   function rebuilds the product from .ahat, .uv and flux.C. Failing both it
%   says so once and draws the lines alone.
%
%   ---------------------------------------------------------------------
%   WHY 45 DAYS
%   ---------------------------------------------------------------------
%   Spring/neap is the M2-S2 beat: 1/(2 - 1.9323) cpd = 14.765 d. One cycle
%   in the window cannot distinguish a fortnightly modulation from a single
%   event that happened to rise and fall; two can be a coincidence of two
%   events. Three cycles -- 44.3 d, so 45 -- is the shortest window in which
%   a regular spacing is a claim rather than an impression, and it is short
%   enough that the figure still resolves individual tidal-band bursts.
%
%   Push it to 60-90 d if the fortnight is there and the question becomes
%   whether its phase drifts; shorten it only to zoom a feature already
%   found in the longer view.
%
%   ---------------------------------------------------------------------
%   THE 14.765-DAY RULER
%   ---------------------------------------------------------------------
%   opt.ruler draws dotted verticals every 14.765 d. Its PHASE IS ARBITRARY
%   -- it is anchored at t0, not at syzygy, because nothing in this
%   repository knows the lunar phase. It is a measuring stick for SPACING
%   only: if the peaks in a panel keep step with the ticks, the spacing is
%   fortnightly; if they drift steadily across them, it is not. Do not read
%   a tick as a predicted spring tide.
%
%   A barotropic flux modulated at the beat is forcing. A baroclinic flux
%   that is not, at a site whose barotropic flux is, means the local
%   conversion is not what is setting the arriving wave.
%
%   author TEJ, with help from Claude
%   12 November 2026

%% ---- arguments -------------------------------------------------------
if nargin < 1, R  = []; end
if nargin < 2, t0 = []; end
if nargin < 3, opt = struct(); end

def = struct('ndays',45, 't_end',[], 'smooth_days',3, 'ruler',true, ...
             'raw',true, 'sameY',true, 'kW',true, 'fignum',21, 'cachedir','');
for f = fieldnames(def)'
    % isempty(false) is false, so a logical option set to false survives this.
    if ~isfield(opt, f{1}) || isempty(opt.(f{1})), opt.(f{1}) = def.(f{1}); end
end
if isempty(opt.cachedir)
    opt.cachedir = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data');
end

if isempty(R), R = {'F2','F3'}; end
if isstruct(R), R = num2cell(R); end

% A cellstr of names means load; a cell of structs is already loaded.
if iscellstr(R) %#ok<ISCLSTR>
    names = R;  R = cell(1, numel(names));
    for i = 1:numel(names)
        fn = fullfile(opt.cachedir, sprintf('flux_%s.mat', names{i}));
        if ~isfile(fn)
            error('pflux_uvp_series:nocache', ...
                ['No %s. Run turner_fleat_IT_flux section 1 for %s first, ' ...
                 'or pass the in-memory R.'], fn, names{i});
        end
        R{i} = load(fn).out;
        fprintf('[%s] from %s\n', names{i}, fn);
    end
end

% ---- keep only the ones carrying a flux ------------------------------
keep = false(1, numel(R));
for i = 1:numel(R)
    r = R{i};
    keep(i) = isstruct(r) && isfield(r,'ok') && r.ok && ...
              isfield(r,'flux') && isstruct(r.flux) && ...
              isfield(r.flux,'Fx') && ~isempty(r.time);
end
if ~any(keep)
    error('pflux_uvp_series:noflux', 'None of the %d inputs carried a flux.', numel(R));
end
if ~all(keep)
    for i = find(~keep)
        who = 'input';
        if isstruct(R{i}) && isfield(R{i},'name'), who = R{i}.name; end
        fprintf(2, '  %s has no flux -- dropped\n', who);
    end
end
R  = R(keep);
nm = numel(R);
h.names = cell(1,nm);
for i = 1:nm
    if isfield(R{i},'name'), h.names{i} = R{i}.name; else, h.names{i} = sprintf('#%d',i); end
end

%% ---- the window ------------------------------------------------------
% Default start is the first hour EVERY mooring has, so the leftmost panel
% edge is not one record's gap.
tfirst = -inf;  tlast = inf;
for i = 1:nm
    t = R{i}.time(:)';
    tfirst = max(tfirst, t(1));  tlast = min(tlast, t(end));
end
if isempty(t0), t0 = tfirst; end
if ~isempty(opt.t_end)
    t1 = opt.t_end;
else
    t1 = t0 + opt.ndays;
end

if t0 < tfirst - 1e-6 || t1 > tlast + 1e-6
    warning('pflux_uvp_series:window', ...
        ['Asked for %s to %s but the moorings overlap only %s to %s. ' ...
         'The panels are padded, not extrapolated.'], ...
        datestr(t0,'dd-mmm-yyyy'), datestr(t1,'dd-mmm-yyyy'), ...
        datestr(tfirst,'dd-mmm-yyyy'), datestr(tlast,'dd-mmm-yyyy'));
end
h.t0 = t0;  h.t1 = t1;
ncyc = (t1 - t0)/14.765;

%% ---- gather the series -----------------------------------------------
sc  = 1;  ulab = 'W m^{-1}';
if opt.kW, sc = 1e-3;  ulab = 'kW m^{-1}'; end

% ser(k,i): panel k (1 = u, 2 = v), mooring i.
ser    = cell(2, nm);
rawsrc = repmat({'none'}, 1, nm);
for i = 1:nm
    t  = R{i}.time(:)';
    in = t >= t0 - 1e-6 & t <= t1 + 1e-6;
    dt = 1;
    if numel(t) > 1, dt = median(diff(t))*24; end

    Fx = R{i}.flux.Fx(:)' * sc;
    Fy = R{i}.flux.Fy(:)' * sc;

    % the instantaneous product, if this record can supply it at all
    [Rx, Ry, src] = raw_uvp(R{i});
    if numel(Rx) ~= numel(Fx) || numel(Ry) ~= numel(Fy)
        Rx = nan(size(Fx));  Ry = nan(size(Fy));  src = 'none';
    else
        Rx = Rx * sc;  Ry = Ry * sc;
    end
    rawsrc{i} = src; %#ok<AGROW>

    ser{1,i} = struct('t', t(in), 'y', Fx(in), 'yr', Rx(in), 'dt', dt);
    ser{2,i} = struct('t', t(in), 'y', Fy(in), 'yr', Ry(in), 'dt', dt);
end

% One note, not one per panel per mooring.
if opt.raw
    miss = h.names(strcmp(rawsrc, 'none'));
    if ~isempty(miss)
        fprintf(2, ['  no instantaneous flux for %s -- flux.Fx_inst is absent. ' ...
                    'Rerun the mooring (compute_tidal_flux saves it as of ' ...
                    '12 Nov 2026) or pass an opt.keep = ''full'' out. ' ...
                    'Drawing the lines only.\n'], strjoin(miss, ', '));
    end
    reb = h.names(strcmp(rawsrc, 'rebuilt'));
    if ~isempty(reb)
        fprintf('  %s: instantaneous flux rebuilt from .ahat, .uv and flux.C\n', ...
                strjoin(reb, ', '));
    end
end
ptitle = {'\langle u''p''\rangle    eastward', ...
          '\langle v''p''\rangle    northward'};
pcomp  = {'u','v'};

% One ruler for both panels: u and v are the same quantity on two axes.
% The dots are included in the limits when they are drawn: clipping them
% would make the scatter look narrower than it is, which is the one thing
% the dots are there to show. opt.raw = false gets the tight scale back.
ylo = 0;  yhi = 0;
for k = 1:2
    for i = 1:nm
        y = ser{k,i}.y(isfinite(ser{k,i}.y));
        if ~isempty(y), ylo = min(ylo, min(y));  yhi = max(yhi, max(y)); end
        if opt.raw
            y = ser{k,i}.yr(isfinite(ser{k,i}.yr));
            if ~isempty(y), ylo = min(ylo, min(y));  yhi = max(yhi, max(y)); end
        end
    end
end
pad = 0.06*max(yhi - ylo, eps);
YL  = [ylo - pad, yhi + pad];

%% ---- draw ------------------------------------------------------------
cM = lines(max(nm,2));

h.fig = figure(opt.fignum); clf
set(h.fig, 'Color','w', 'Position',[80 120 980 580])
T = tiledlayout(h.fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

h.ax = gobjects(1,2);
hl1  = gobjects(1,0);  lab1 = cell(1,0);
for k = 1:2
    ax = nexttile(T);  h.ax(k) = ax;  hold(ax,'on');  box(ax,'on')

    % the ruler first, so the data sits on top of it
    if opt.ruler
        for tk = t0 : 14.765 : t1
            xline(ax, tk, ':', 'Color',[.6 .6 .6], 'LineWidth',0.8);
        end
    end
    yline(ax, 0, 'k-', 'LineWidth',0.6);

    % The legend names EVERY line, not one per mooring: the three traces are
    % three different processings of the same series and which is which is the
    % one thing a reader has to know before reading anything off the panel.
    % That belongs next to the line it describes rather than in a subtitle the
    % eye has already left. It runs to three entries per mooring and it is
    % clunky, which is the price of the figure being self-describing.
    hl  = gobjects(1,0);
    lab = cell(1,0);
    any_data = false;
    for i = 1:nm
        s = ser{k,i};
        if isempty(s.t) || ~any(isfinite(s.y))
            hl(end+1) = plot(ax, NaN, NaN, '-', 'Color',cM(i,:), 'LineWidth',1.7); %#ok<AGROW>
            lab{end+1} = sprintf('%s  (no data)', h.names{i}); %#ok<AGROW>
            continue
        end
        any_data = true;

        % Dots first and faint: they are the widest thing on the panel and
        % would bury the lines they exist to be compared against.
        if opt.raw && any(isfinite(s.yr))
            plot(ax, s.t, s.yr, '.', 'Color',[cM(i,:) 0.18], 'MarkerSize',3.5);
            % A 3.5-pt dot at 18% alpha is invisible in a legend key, so the
            % entry gets an opaque proxy at the marker size a key needs. Same
            % colour, same marker, nothing on the axes drawn twice.
            hl(end+1)  = plot(ax, NaN, NaN, '.', 'Color',cM(i,:), 'MarkerSize',12); %#ok<AGROW>
            lab{end+1} = sprintf('%s  instantaneous (2\\omega, not a flux)', ...
                                 h.names{i}); %#ok<AGROW>
        end

        % The thin line is the 25-h phase average as computed; the heavy line
        % is the running mean. Same hue, so colour never stops meaning mooring.
        thin = plot(ax, s.t, s.y, '-', 'Color',[cM(i,:) 0.40], 'LineWidth',0.9);
        % The boxcar compute_tidal_flux actually used, in hours, read back off
        % cfg rather than hard-coded: 2 periods at hourly sampling is 25 h, but
        % cfg.flux_periods is a knob and a label that lies is worse than none.
        % An empty flux_periods means the function's own default of 2.
        nper = 25;
        if isfield(R{i},'cfg') && all(isfield(R{i}.cfg, {'M2_cpd','dt_hr'}))
            per = 2;
            if isfield(R{i}.cfg,'flux_periods') && ~isempty(R{i}.cfg.flux_periods)
                per = R{i}.cfg.flux_periods;
            end
            nper = round(per * (24/R{i}.cfg.M2_cpd) / R{i}.cfg.dt_hr) * R{i}.cfg.dt_hr;
        end
        tlab = sprintf('%s  %g-h phase average', h.names{i}, nper);
        cov  = nnz(isfinite(s.y))/numel(s.y);
        if cov < 0.995, tlab = sprintf('%s  (%.0f%% cover)', tlab, 100*cov); end

        if opt.smooth_days > 0
            nw   = max(3, round(opt.smooth_days*24/s.dt));
            heav = plot(ax, s.t, movmean(s.y, nw, 'omitnan'), '-', ...
                        'Color',cM(i,:), 'LineWidth',1.8);
            hl(end+1)  = thin;  lab{end+1} = tlab; %#ok<AGROW>
            hl(end+1)  = heav; %#ok<AGROW>
            lab{end+1} = sprintf('%s  %g-d running mean', ...
                                 h.names{i}, opt.smooth_days); %#ok<AGROW>
        else
            set(thin, 'Color', cM(i,:), 'LineWidth', 1.3);
            hl(end+1) = thin;  lab{end+1} = tlab; %#ok<AGROW>
        end
    end

    xlim(ax, [t0 t1])
    if opt.sameY && diff(YL) > 0, ylim(ax, YL); end
    grid(ax,'on')
    ylabel(ax, ulab)
    title(ax, ptitle{k}, 'FontWeight','normal', 'HorizontalAlignment','left', ...
          'Units','normalized', 'Position',[0.012 0.98 0])
    % ONE legend for the figure, from the first panel, docked under the stack.
    % Three entries per mooring will not fit inside an axes without covering
    % the series, and the labels are identical in both panels -- u and v come
    % off the same fit, so even the coverage figure matches. Drawing it twice
    % would be two thirds of a tile of repetition.
    if k == 1, hl1 = hl;  lab1 = lab; end
    if ~any_data
        text(ax, 0.5, 0.45, 'no data in this window', 'Units','normalized', ...
             'HorizontalAlignment','center', 'Color',[.5 .5 .5])
    end

    if k < 2
        set(ax, 'XTickLabel', []);
    else
        if t1 - t0 <= 75
            datetick(ax, 'x', 'dd mmm', 'keeplimits');
        else
            datetick(ax, 'x', 'mmm yyyy', 'keeplimits');
        end
        xlabel(ax, sprintf('%s to %s', datestr(t0,'dd mmm yyyy'), ...
                                       datestr(t1,'dd mmm yyyy')))
    end
end

if ~isempty(hl1)
    lg = legend(h.ax(1), hl1, lab1, 'Box','off', 'FontSize',8.5, ...
                'NumColumns', max(nm,1));
    % Docking to the layout needs R2020b. Older MATLAB keeps it on the axes,
    % outside them, which costs height but never covers the data.
    try
        lg.Layout.Tile = 'south';
    catch
        set(lg, 'Location', 'southoutside');
    end
    h.lgd = lg;
end

% What each trace IS now lives in the legend, beside the line it describes.
% The subtitle keeps only what is true of the whole figure.
sub = sprintf('%s  |  %.0f d = %.1f spring/neap cycles', ...
              strjoin(h.names, ', '), t1-t0, ncyc);
if opt.ruler
    %sub = sprintf('%s  |  dotted: 14.765 d, phase arbitrary', sub);
end
title(T, {'depth-integrated baroclinic M_2 energy flux', sub}, 'FontWeight','bold')

%% ---- what the console should say -------------------------------------
fprintf('\nwindow %s to %s  (%.0f d, %.1f spring/neap cycles)\n', ...
        datestr(t0,'dd-mmm-yyyy'), datestr(t1,'dd-mmm-yyyy'), t1-t0, ncyc);
for k = 1:2
    for i = 1:nm
        s = ser{k,i};
        nmlab = sprintf('%s %s', h.names{i}, pcomp{k});
        if isempty(s.t), fprintf('  %-6s no data in this window\n', nmlab); continue, end
        g = s.y(isfinite(s.y));
        if isempty(g), fprintf('  %-6s all NaN in this window\n', nmlab); continue, end
        fprintf('  %-6s mean %+7.3f   rms %6.3f   range %+7.3f to %+7.3f  %s\n', ...
                nmlab, mean(g), sqrt(mean(g.^2)), min(g), max(g), ulab);
    end
end
fprintf(['\nRead the heavy lines: peaks a fortnight apart, keeping step with the\n' ...
         'dotted ticks, are spring/neap. Peaks that drift across them are not.\n\n']);
end

% =======================================================================
function [Rx, Ry, src] = raw_uvp(r)
%RAW_UVP  Instantaneous int u'p' dz [W/m], by whichever route the record allows.
%
%   src  'saved'    straight off flux.Fx_inst / .Fy_inst
%        'rebuilt'  from .ahat.band, .uv.uhat_M2/.vhat_M2 and flux.C
%        'none'     neither available; Rx and Ry come back empty
%
%   The rebuild is the same single sum compute_tidal_flux does,
%
%       int u'p' dz = sum_i uhat_i(t) ahat_i(t) C_i
%
%   with the +1 on uhat because ITS ROW 1 IS BAROTROPIC while ahat and C are
%   baroclinic-only. It is a fallback for an out that predates Fx_inst, not a
%   second implementation to maintain: if the two ever disagree, the saved one
%   is right and this one has drifted.

Rx = [];  Ry = [];  src = 'none';
if ~isfield(r,'flux') || ~isstruct(r.flux), return, end

if isfield(r.flux,'Fx_inst') && isfield(r.flux,'Fy_inst') && ...
        ~isempty(r.flux.Fx_inst)
    Rx = r.flux.Fx_inst(:)';  Ry = r.flux.Fy_inst(:)';  src = 'saved';
    return
end

% ---- the keep = 'full' route ----------------------------------------
if ~isfield(r,'ahat') || ~isfield(r,'uv') || ~isfield(r.flux,'C'), return, end
if ~isfield(r.ahat,'band') || ~isfield(r.uv,'uhat_M2') || ~isfield(r.uv,'vhat_M2')
    return
end

a = r.ahat.band;  u = r.uv.uhat_M2;  v = r.uv.vhat_M2;  C = r.flux.C(:);
nf = min([numel(C), size(a,1), size(u,1)-1, size(v,1)-1]);
if nf < 1 || size(a,2) ~= numel(r.time), return, end

Rx = zeros(1, size(a,2));  Ry = Rx;
for q = 1:nf
    Rx = Rx + u(q+1,:) .* a(q,:) * C(q);
    Ry = Ry + v(q+1,:) .* a(q,:) * C(q);
end
src = 'rebuilt';
end
