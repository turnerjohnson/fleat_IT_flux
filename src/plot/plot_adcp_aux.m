function out = plot_adcp_aux(m, opt)
% PLOT_ADCP_AUX  The four non-velocity ADCP records, and what they really are.
%
%   out = plot_adcp_aux(m)
%   out = plot_adcp_aux(m, opt)
%
%   m    mooring number 1-5. Default 1.
%   opt .datadir  where FLEAT_M_ADCP_8m_BindDp.mat lives
%       .tidx     time indices to plot. [] (default) for the full record
%       .target   decimate to about this many columns. Default 4000
%
%   out .time .depth .Tup .Tdn   the plotted series
%       .offset .step .spread    the clock, in minutes
%       .stats                   the table view, printed as well as returned
%       .fig                     the figure handle
%
%   THE FOUR VARIABLES, NAMED HONESTLY
%   FLEAT_M_ADCP_8m_BindDp.mat carries four [1 x nt] records per mooring
%   besides U and V, and two of the four are misnamed in the file:
%
%     dnADCPdepth   REAL. Depth of the down-looking head, m. This is the
%                   knockdown record -- the mooring leaning over in the flow.
%                   Median 57-89 m depending on mooring, with excursions to
%                   186 m at M2. Everything hanging below it moves with it.
%
%     upADCPtilt    NOT TILT. Instrument TEMPERATURE in degC, from the
%     dnADCPtilt    up- and down-looking heads. build_dod_mooring establishes
%                   this: r = +0.98 to +0.99 against TZdepth at the
%                   instrument's own depth, bias < 0.33 degC, and the values
%                   run 11-31 degC, which is a temperature and not a degree of
%                   tilt. They are stored as F.adcp.Tup and F.adcp.Tdn for
%                   exactly that reason. Do NOT use either as a knockdown or
%                   tilt flag; dnADCPdepth is the knockdown record.
%
%     ADCPtm        the sample times, datenum. Plotted here as MINUTES PAST
%                   THE HOUR, because that is the only thing in it that varies
%                   -- and it does not vary either: it is a flat line at 14 to
%                   55 min depending on the mooring. See
%                   PLOT_ADCP_TIME_MISALIGN for the five-mooring version of
%                   that fact, and WHICH WAY TO MOVE THE ADCP in
%                   build_dod_mooring for why it matters.
%
%   WHY Tup AND Tdn SHARE A Y AXIS
%   They are the same quantity measured a few metres apart on the same wire,
%   so the comparison between them is the content: at M1, M2, M3 and M5 they
%   track each other within ~0.15 degC of median, and at M4 they do not --
%   Tup sits 3.6 degC colder in the median AND is 56% NaN, the same upper-bin
%   dropout that build_dod_mooring warns about at M4. Autoscaling each panel
%   separately would hide that by construction, so both panels are forced onto
%   one set of limits.
%
%   The depth panel is NOT on those limits, and is drawn reversed (positive
%   down) like every other depth axis in this directory.

%% ---- arguments ------------------------------------------------------
if nargin < 1 || isempty(m), m = 1; end
if nargin < 2 || ~isstruct(opt), opt = struct(); end
if ~isfield(opt,'datadir') || isempty(opt.datadir)
    opt.datadir = '/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/';
end
if ~isfield(opt,'tidx'),                            opt.tidx   = [];   end
if ~isfield(opt,'target') || isempty(opt.target),   opt.target = 4000; end
if ~isnumeric(m) || ~isscalar(m) || m < 1 || m > 5
    error('plot_adcp_aux:mooring', 'm is a single mooring number 1-5. Got %s.', mat2str(m));
end

fa = fullfile(opt.datadir, 'FLEAT_M_ADCP_8m_BindDp.mat');
if ~isfile(fa), error('plot_adcp_aux:nofile', 'Missing: %s', fa); end

C.surface = [252 252 251]/255;
C.ink     = [ 11  11  11]/255;
C.ink2    = [ 82  81  78]/255;
C.muted   = [154 164 173]/255;
C.series  = [ 42 120 214;             % 1 blue     depth
              235 104  52;            % 2 orange   Tup
               27 175 122;            % 3 aqua     Tdn
              237 161   0]/255;       % 4 yellow   clock

%% ---- load -----------------------------------------------------------
A = load(fa, sprintf('ADCPtm%d', m), sprintf('dnADCPdepth%d', m), ...
             sprintf('upADCPtilt%d', m), sprintf('dnADCPtilt%d', m));
t   = double(A.(sprintf('ADCPtm%d',      m))(:)).';
dep = double(A.(sprintf('dnADCPdepth%d', m))(:)).';
Tup = double(A.(sprintf('upADCPtilt%d',  m))(:)).';   % temperature, not tilt
Tdn = double(A.(sprintf('dnADCPtilt%d',  m))(:)).';

off = mod((t - floor(t))*1440, 60);
out.offset = off(1);
out.spread = max(off) - min(off);
out.step   = median(diff(t))*1440;

%% ---- subset and decimate -------------------------------------------
n = numel(t);
if ~isempty(opt.tidx)
    k = opt.tidx(opt.tidx >= 1 & opt.tidx <= n);
    if isempty(k), error('plot_adcp_aux:tidx', 'tidx is outside 1:%d.', n); end
else
    k = 1:n;
end
dec = max(1, floor(numel(k)/opt.target));
if dec > 1
    fprintf('plot_adcp_aux: decimating %d -> %d columns (every %d).\n', ...
        numel(k), numel(1:dec:numel(k)), dec);
    k = k(1:dec:numel(k));
end
t = t(k);  dep = dep(k);  Tup = Tup(k);  Tdn = Tdn(k);  off = off(k);
dt = datetime(t, 'ConvertFrom', 'datenum');

%% ---- figure ---------------------------------------------------------
f  = figure('Color', C.surface, 'Position', [100 80 1000 780]);
tl = tiledlayout(f, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% the two temperatures share limits -- see WHY Tup AND Tdn SHARE A Y AXIS
Tall = [Tup Tdn];  Tall = Tall(isfinite(Tall));
% headroom at the top of both panels so the in-panel labels never sit on data
Tlim = [min(Tall) - 0.03*range(Tall), max(Tall) + 0.22*range(Tall)];
dlim = [min(dep) - 0.22*range(dep),   max(dep) + 0.03*range(dep)];

ser = { dep, 'dnADCPdepth', 'head depth (m)',       1, true,  dlim ; ...
        Tup, 'upADCPtilt',  'T_{up} (\circC)',      2, false, Tlim ; ...
        Tdn, 'dnADCPtilt',  'T_{dn} (\circC)',      3, false, Tlim ; ...
        off, 'ADCPtm',      'min past hour',        4, false, [] };

ax = gobjects(1, 4);
for i = 1:4
    ax(i) = nexttile(tl);
    hold(ax(i), 'on');
    set(ax(i), 'Color', C.surface, 'XColor', C.ink2, 'YColor', C.ink2, ...
        'FontSize', 9, 'TickDir', 'out', 'Box', 'off');
    y = ser{i,1};
    plot(ax(i), dt, y, '-', 'Color', C.series(ser{i,4},:), 'LineWidth', 0.8);
    ylabel(ax(i), ser{i,3}, 'Color', C.ink, 'FontSize', 10);
    if ~isempty(ser{i,6}), ylim(ax(i), ser{i,6}); end
    if ser{i,5}, set(ax(i), 'YDir', 'reverse'); end
    xlim(ax(i), [dt(1) dt(end)]);
    if i < 4, set(ax(i), 'XTickLabel', []); end

    % name the field, and say what it is where the file's name is wrong
    lab = ser{i,2};
    if i == 2 || i == 3
        lab = sprintf('%s  --  TEMPERATURE, not tilt', ser{i,2});
    elseif i == 4
        lab = sprintf('%s  --  %.2f min past the hour, step %.3f min, spread %.3f min', ...
            ser{i,2}, out.offset, out.step, out.spread);
    end
    % panel 4's line sits wherever the offset is, so its label goes to the
    % floor of the panel instead of the ceiling
    if i == 4
        text(ax(i), 0.006, 0.07, lab, 'Units', 'normalized', 'Color', C.ink, ...
            'FontSize', 9, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
    else
        text(ax(i), 0.006, 0.93, lab, 'Units', 'normalized', 'Color', C.ink, ...
            'FontSize', 9, 'FontWeight', 'bold', 'VerticalAlignment', 'top');
    end

    nnan = nnz(~isfinite(y));
    if nnan > 0
        text(ax(i), 0.994, 0.93, sprintf('%.0f%% NaN', 100*nnan/numel(y)), ...
            'Units', 'normalized', 'Color', C.ink2, 'FontSize', 8.5, ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
    end
end
linkaxes(ax, 'x');
ylim(ax(4), [0 60]);  yticks(ax(4), 0:15:60);
xlabel(ax(4), 'time (UTC)', 'Color', C.ink);
ax(4).XAxis.SecondaryLabel.Visible = 'off';

title(tl, sprintf('M%d ADCP: the four records that are not velocity', m), ...
    'Color', C.ink, 'FontSize', 13, 'FontWeight', 'bold');
subtitle(tl, ['dnADCPdepth is the knockdown; the two ''tilt'' fields are ' ...
    'instrument temperature; ADCPtm is the clock'], ...
    'Color', C.ink2, 'FontSize', 10);

out.time = t;  out.depth = dep;  out.Tup = Tup;  out.Tdn = Tdn;
out.fig  = f;

%% ---- the table view -------------------------------------------------
fprintf('\nplot_adcp_aux: M%d, %d of %d samples plotted\n', m, numel(t), n);
fprintf('  field          min      med      max     NaN%%\n');
rows = {'dnADCPdepth', dep; 'upADCPtilt', Tup; 'dnADCPtilt', Tdn};
out.stats = nan(3,4);
for i = 1:3
    y = rows{i,2};
    out.stats(i,:) = [min(y,[],'omitnan') median(y,'omitnan') max(y,[],'omitnan') ...
                      100*nnz(~isfinite(y))/numel(y)];
    fprintf('  %-12s %8.2f %8.2f %8.2f %8.1f\n', rows{i,1}, out.stats(i,:));
end
fprintf('  ADCPtm       %.2f min past the hour, step %.3f min, spread %.3f min over the record\n', ...
    out.offset, out.step, out.spread);
dT = median(Tup - Tdn, 'omitnan');
fprintf('  median Tup - Tdn = %+.2f degC', dT);
if abs(dT) > 1
    fprintf('  <-- the two heads disagree; at M4 this is the 56%% upper-bin dropout\n');
else
    fprintf('  (same water, a few metres apart)\n');
end
end
