function out = plot_adcp_time_misalign(opt)
% PLOT_ADCP_TIME_MISALIGN  When each DOD ADCP actually samples, M1-M5.
%
%   out = plot_adcp_time_misalign()
%   out = plot_adcp_time_misalign(opt)
%
%   opt .datadir  where FLEAT_M_ADCP_8m_BindDp.mat lives
%       .hours    length of the window, h. Default 8
%       .t0       datenum of the window start. Default [] -- one day into the
%                 common overlap, which is arbitrary and does not matter: the
%                 offsets are the same everywhere in the record
%       .ax       axes to draw into. Default [] -- a new figure
%
%   out .t       [5 x 1] cell, the sample times in the window (datenum)
%       .offset  [1 x 5] minutes past the hour, per mooring
%       .spread  [1 x 5] max-min of that offset over the WHOLE record, min
%       .step    [1 x 5] median sample interval, minutes
%       .ax      the axes handle
%
%   WHAT THIS DRAWS, AND WHY IT IS ONLY MARKERS
%   One row per mooring, one marker per sample, against the exact UTC hours as
%   dotted verticals. No velocity: this figure is about the clock, and a u
%   trace invites the eye to compare shapes between moorings when the only
%   claim being made is WHEN each instrument wrote a sample. The markers land
%   where they land, 14 to 55 minutes past the hour, and never on one.
%
%   The offsets are read from ADCPtm, not assumed, and they are CONSTANT --
%   out.spread is the max-min over all ~7540 samples and it is 0.00 min at four
%   moorings, 0.02 min at M4. So there is no drift to model: one fixed shift
%   per mooring is the whole story, which is what lets build_dod_mooring put u
%   and v on an exact-hour clock with a single interpolation (see WHICH WAY TO
%   MOVE THE ADCP in its header -- the choice there costs either ~3% of M2
%   amplitude or up to 14.5 deg of M2 phase, and this is the picture of why
%   that choice exists at all).
%
%   M1 and M5 are both ~35 min past the hour, so their markers sit on the same
%   verticals. That is the data, not overplotting. The worst pair is M3 (14
%   min) against M4 (55 min), 41 minutes apart.
%
%   Colour is backed by marker shape and by the y-axis labels, so identity
%   never rests on hue alone -- which it cannot here, since three of these five
%   validated categorical steps sit under 3:1 contrast on a light surface.

%% ---- arguments ------------------------------------------------------
if nargin < 1 || ~isstruct(opt), opt = struct(); end
if ~isfield(opt,'datadir') || isempty(opt.datadir)
    opt.datadir = '/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/';
end
if ~isfield(opt,'hours') || isempty(opt.hours), opt.hours = 8;  end
if ~isfield(opt,'t0'),                          opt.t0    = []; end
if ~isfield(opt,'ax'),                          opt.ax    = []; end

fa = fullfile(opt.datadir, 'FLEAT_M_ADCP_8m_BindDp.mat');
if ~isfile(fa)
    error('plot_adcp_time_misalign:nofile', 'Missing: %s', fa);
end

% palette -- validated categorical slots 1-5, light surface
C.surface = [252 252 251]/255;
C.ink     = [ 11  11  11]/255;
C.ink2    = [ 82  81  78]/255;
C.muted   = [154 164 173]/255;
C.series  = [ 42 120 214;             % 1 blue     M1
              235 104  52;            % 2 orange   M2
               27 175 122;            % 3 aqua     M3
              237 161   0;            % 4 yellow   M4
              232 123 164]/255;       % 5 magenta  M5
mk = {'o','s','^','d','v'};           % secondary encoding, not decoration

%% ---- the clocks, read rather than assumed ---------------------------
nm = 5;
A  = load(fa, 'ADCPtm1','ADCPtm2','ADCPtm3','ADCPtm4','ADCPtm5');
ta = cell(1, nm);
out.offset = nan(1, nm);  out.spread = nan(1, nm);  out.step = nan(1, nm);
for m = 1:nm
    t   = double(A.(sprintf('ADCPtm%d', m))(:)).';
    off = mod((t - floor(t))*1440, 60);
    ta{m} = t;
    out.offset(m) = off(1);
    out.spread(m) = max(off) - min(off);
    out.step(m)   = median(diff(t))*1440;
end

%% ---- the window -----------------------------------------------------
span = opt.hours/24;
if isempty(opt.t0)
    t0 = max(cellfun(@(t) t(1), ta)) + 1;
else
    t0 = opt.t0;
end
t1  = t0 + span;
hrs = ceil(t0*24)/24 : 1/24 : t1;                 % the exact UTC hours
xl  = datetime([t0 t1], 'ConvertFrom', 'datenum');

%% ---- the axes -------------------------------------------------------
if isempty(opt.ax)
    f  = figure('Color', C.surface, 'Position', [100 300 980 300]);
    ax = axes(f);                                                        %#ok<LAXES>
else
    ax = opt.ax;
end
hold(ax, 'on');
set(ax, 'Color', C.surface, 'XColor', C.ink2, 'YColor', C.ink2, ...
    'FontSize', 10, 'TickDir', 'out', 'Box', 'off');

for h = datetime(hrs, 'ConvertFrom', 'datenum')
    xline(ax, h, ':', 'Color', C.muted, 'LineWidth', 0.6, 'HandleVisibility', 'off');
end

for m = 1:nm
    g  = ta{m} >= t0 & ta{m} < t1;
    dt = datetime(ta{m}(g), 'ConvertFrom', 'datenum');
    out.t{m,1} = ta{m}(g);
    plot(ax, dt, (nm + 1 - m)*ones(size(dt)), 'LineStyle', 'none', ...
        'Marker', mk{m}, 'MarkerSize', 8, 'LineWidth', 1.2, ...
        'MarkerFaceColor', C.series(m,:), 'MarkerEdgeColor', C.surface, ...
        'DisplayName', sprintf('M%d  (+%.0f min)', m, out.offset(m)));
end

xlim(ax, xl);  ylim(ax, [0.4 nm + 0.6]);
yticks(ax, 1:nm);
yticklabels(ax, arrayfun(@(m) sprintf('M%d', m), nm:-1:1, 'UniformOutput', false));
xlabel(ax, sprintf('time, %s UTC   (dotted lines = exact UTC hours)', ...
    datestr(t0, 'yyyy-mm-dd')), 'Color', C.ink);
ax.XAxis.SecondaryLabel.Visible = 'off';
legend(ax, 'Location', 'eastoutside', 'TextColor', C.ink, ...
    'Color', C.surface, 'EdgeColor', C.ink2, 'FontSize', 9);

out.ax = ax;
out.t0 = t0;

%% ---- the table view -------------------------------------------------
fprintf('\nplot_adcp_time_misalign: %g h from %s UTC\n', ...
    opt.hours, datestr(t0, 'yyyy-mm-dd HH:MM'));
fprintf('  moor  offset(min)  spread(min)  step(min)  n in window\n');
for m = 1:nm
    fprintf('   M%d   %9.2f  %11.3f  %9.3f  %11d\n', m, out.offset(m), ...
        out.spread(m), out.step(m), numel(out.t{m}));
end
[~, a] = min(out.offset);  [~, b] = max(out.offset);
dmin = out.offset(b) - out.offset(a);
fprintf(['  worst pair M%d vs M%d: %.1f min apart = %.1f deg of M2 (12.42 h) ' ...
         'if u,v are snapped to the nearest hour\n'], ...
    a, b, dmin, 360*dmin/(12.42*60));
end
