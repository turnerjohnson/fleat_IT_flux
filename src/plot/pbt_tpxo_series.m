function h = pbt_tpxo_series(R, t0, opt)
%PBT_TPXO_SERIES  Barotropic M2 current through time, same window as the flux.
%
%   h = pbt_tpxo_series()                        % F2 and F3
%   h = pbt_tpxo_series({'F2','F3'}, t0)
%   h = pbt_tpxo_series(R, t0, opt)              % R from turner_fleat_IT_flux
%
%   The companion to pflux_uvp_series, drawn the same way -- u on top, v
%   below, colour carrying the mooring, one ruler, one window -- so the two
%   figures can be read against each other on screen.
%
%   R    cellstr of mooring names (their configs supply lat, lon and H), OR
%        cell/struct array of fleat_IT_flux outputs, OR [] for {'F2','F3'}.
%   t0   window start, datenum. [] means 01-Sep-2016, a date inside both
%        records; nothing here needs the mooring record at all.
%   opt .ndays   window length [d]. Default 45, to match the flux figure.
%       .t_end   datenum, overrides ndays.
%       .dt_hr   prediction step [h]. Default 0.25 -- see WHY 15 MINUTES.
%       .H       scalar or per-mooring vector, overriding each cfg.H [m].
%       .file    the atlas .nc, passed to tpxo_m2_bt.
%       .ruler   true (default) marks 14.765-d intervals from t0.
%       .cms     true (default) plots cm/s rather than m/s.
%       .fignum  default 22.
%
%   h .fig .ax(2) .lgd .names .B {1 x nm} the tpxo_m2_bt outputs .t0 .t1
%
%   ---------------------------------------------------------------------
%   READ THE HEADER OF TPXO_M2_BT BEFORE READING THIS FIGURE
%   ---------------------------------------------------------------------
%   Two things about it are not visible on the axes:
%
%   1. The atlas holds TRANSPORT. The velocity is transport/H, and H here is
%      the mooring's recorded depth, not TPXO's own h at that node, because
%      grid_tpxo10_atlas_30_v2.nc is not in data/. That is a pure SCALE error
%      on the amplitude -- the shape and phase are unaffected -- so the figure
%      cannot show it. The subtitle prints the H it used.
%
%   2. It is M2 ALONE, so THERE IS NO SPRING/NEAP HERE. The fortnight is the
%      M2-S2 beat and S2 is in a different file. The dashed envelope is drawn
%      flat on purpose: it is the figure saying so. Anything fortnightly in
%      the baroclinic flux cannot be compared against this curve's envelope,
%      only against its phase. Fetch u_s2_tpxo10_atlas_30_v2.nc and sum the
%      two predictions to get a barotropic envelope worth comparing.
%
%   ---------------------------------------------------------------------
%   WHY 15 MINUTES, AND WHAT 45 DAYS LOOKS LIKE
%   ---------------------------------------------------------------------
%   M2 is 12.42 h, so a 45-day window holds 87 cycles. At any linewidth that
%   renders as a solid band, which is honest -- the envelope IS flat -- but it
%   is not a waveform. Use opt.ndays = 3 to 5 to see the oscillation itself,
%   and keep 45 for the envelope comparison against the flux figure.
%
%   The prediction is on its own 15-minute clock rather than the record's
%   hourly one: hourly sampling of a 12.42-h signal is only 12 points per
%   cycle and the drawn peaks would wander by a few percent with the aliasing
%   rather than sitting on the true envelope.
%
% author TEJ, with help from Claude
% 12 November 2026

%% ---- arguments -------------------------------------------------------
if nargin < 1, R  = []; end
if nargin < 2, t0 = []; end
if nargin < 3, opt = struct(); end

def = struct('ndays',45, 't_end',[], 'dt_hr',0.25, 'H',[], 'file','', ...
             'ruler',true, 'cms',true, 'fignum',22);
for f = fieldnames(def)'
    if ~isfield(opt, f{1}) || isempty(opt.(f{1})), opt.(f{1}) = def.(f{1}); end
end

if isempty(R), R = {'F2','F3'}; end
if isstruct(R), R = num2cell(R); end
if isempty(t0), t0 = datenum(2016,9,1); end
if ~isempty(opt.t_end), t1 = opt.t_end; else, t1 = t0 + opt.ndays; end

% ---- site list: names go through the configs, outs carry their own -----
nm = numel(R);
h.names = cell(1,nm);  lonq = nan(1,nm);  latq = nan(1,nm);  Hq = nan(1,nm);
for i = 1:nm
    if ischar(R{i}) || isstring(R{i})
        nmi = char(R{i});
        cfg = feval([lower(nmi) '_config']);
        h.names{i} = nmi;  lonq(i) = cfg.lon;  latq(i) = cfg.lat;  Hq(i) = cfg.H;
    elseif isstruct(R{i})
        r = R{i};
        if isfield(r,'name'), h.names{i} = r.name; else, h.names{i} = sprintf('#%d',i); end
        % A failed out still carries its name but has NaN for the position, so
        % fall back to the config rather than dropping the mooring.
        if isfield(r,'lon') && isfinite(r.lon)
            lonq(i) = r.lon;  latq(i) = r.lat;  Hq(i) = r.H;
        else
            cfg = feval([lower(h.names{i}) '_config']);
            lonq(i) = cfg.lon;  latq(i) = cfg.lat;  Hq(i) = cfg.H;
        end
    else
        error('pbt_tpxo_series:input', 'Input %d is neither a name nor a struct.', i);
    end
end
if ~isempty(opt.H)
    if isscalar(opt.H), Hq(:) = opt.H; else, Hq = opt.H(:)'; end
end

%% ---- predict ---------------------------------------------------------
tv  = t0 : opt.dt_hr/24 : t1;
sc  = 1;  ulab = 'm s^{-1}';
if opt.cms, sc = 100;  ulab = 'cm s^{-1}'; end

topt = struct('file', opt.file);
h.B  = cell(1,nm);
for i = 1:nm
    h.B{i} = tpxo_m2_bt(lonq(i), latq(i), tv, Hq(i), topt);
end
h.t0 = t0;  h.t1 = t1;

%% ---- draw ------------------------------------------------------------
cM = lines(max(nm,2));
comp   = {'u','v'};
ptitle = {'u_{BT}    eastward', 'v_{BT}    northward'};

ylim_all = 0;
for i = 1:nm
    ylim_all = max([ylim_all, sc*h.B{i}.amp_u, sc*h.B{i}.amp_v]);
end
YL = [-1.12 1.12] * max(ylim_all, eps);

h.fig = figure(opt.fignum); clf
set(h.fig, 'Color','w', 'Position',[80 120 980 580])
T = tiledlayout(h.fig, 2, 1, 'TileSpacing','compact', 'Padding','compact');

h.ax = gobjects(1,2);
hl1  = gobjects(1,0);  lab1 = cell(1,0);
for k = 1:2
    ax = nexttile(T);  h.ax(k) = ax;  hold(ax,'on');  box(ax,'on')

    if opt.ruler
        for tk = t0 : 14.765 : t1
            xline(ax, tk, ':', 'Color',[.6 .6 .6], 'LineWidth',0.8);
        end
    end
    yline(ax, 0, 'k-', 'LineWidth',0.6);

    hl = gobjects(1,0);  lab = cell(1,0);
    for i = 1:nm
        B = h.B{i};
        if k == 1, y = sc*B.u; else, y = sc*B.v; end
        % The envelope drawn is the one actually ON SCREEN: max|y| carries the
        % nodal factor t_predic applied, which the stored amplitude does not.
        a = max(abs(y));

        hs = plot(ax, tv, y, '-', 'Color',[cM(i,:) 0.75], 'LineWidth',0.6);

        % Flat BY CONSTRUCTION, and that is the point: one constituent has no
        % fortnightly envelope to find. Drawn as two lines rather than a band
        % so it cannot be mistaken for an uncertainty.
        he = plot(ax, [t0 t1], [a a], '--', 'Color',cM(i,:), 'LineWidth',1.2);
        plot(ax, [t0 t1], [-a -a], '--', 'Color',cM(i,:), 'LineWidth',1.2);

        hl(end+1)  = hs; %#ok<AGROW>
        lab{end+1} = sprintf('%s  M_2 barotropic, TPXO10 / H = %.0f m', ...
                             h.names{i}, B.H); %#ok<AGROW>
        hl(end+1)  = he; %#ok<AGROW>
        lab{end+1} = sprintf('%s  envelope \\pm%.2f %s (flat: M_2 alone)', ...
                             h.names{i}, a, ulab); %#ok<AGROW>
    end

    xlim(ax, [t0 t1]);  ylim(ax, YL)
    grid(ax,'on');  ylabel(ax, ulab)
    title(ax, ptitle{k}, 'FontWeight','normal', 'HorizontalAlignment','left', ...
          'Units','normalized', 'Position',[0.012 0.98 0])
    if k == 1, hl1 = hl;  lab1 = lab; end

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
    try
        lg.Layout.Tile = 'south';
    catch
        set(lg, 'Location','southoutside');
    end
    h.lgd = lg;
end

sub = sprintf('%s  |  %.0f d  |  transport/H, H = recorded depth  |  M_2 ONLY: no spring/neap here', ...
              strjoin(h.names, ', '), t1-t0);
title(T, {'barotropic M_2 current, TPXO10\_atlas\_v2', sub}, 'FontWeight','bold')

%% ---- console ---------------------------------------------------------
fprintf('\nTPXO M2 barotropic, %s to %s\n', ...
        datestr(t0,'dd-mmm-yyyy'), datestr(t1,'dd-mmm-yyyy'));
for i = 1:nm
    B = h.B{i};
    fprintf(['  %-3s  node %.4f, %.4f (%.1f km off)  H %.0f m\n' ...
             '       |U| %7.1f m^2/s -> u %5.2f cm/s, G %3.0f deg\n' ...
             '       |V| %7.1f m^2/s -> v %5.2f cm/s, G %3.0f deg\n'], ...
            h.names{i}, B.lon_u, B.lat_u, B.moved_km, B.H, ...
            B.U*1e-4, 100*B.amp_u, B.pha_u, ...
            B.V*1e-4, 100*B.amp_v, B.pha_v);
end
fprintf(['\nThe envelope is flat because this is one constituent. To compare a\n' ...
         'fortnightly envelope against the baroclinic flux, add S2.\n\n']);
end
