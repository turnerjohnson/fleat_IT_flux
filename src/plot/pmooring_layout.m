function h = pmooring_layout(opt)
% PMOORING_LAYOUT  The FLEAT array in three dimensions: who sat where.
%
%   h = pmooring_layout()                       % as-deployed, default box
%   h = pmooring_layout(struct('mode','gridded'))
%   h = pmooring_layout(opt)
%
% Draws all five FLEAT moorings (F1 F2 F3 F4 F6 -- there is no F5) as vertical
% strings inside a lon/lat/depth box, over the Palau bathymetry. The figure is
% a normal MATLAB 3-D axes with rotate3d on, so it is meant to be spun.
%
% TWO MODES, AND THE DIFFERENCE BETWEEN THEM IS THE POINT
%
%   'deployed'  (default) every instrument that went in the water, colour-coded
%               by type. This is the mooring as it was built.
%
%   'gridded'   the levels that actually exist in F*_fullfield.mat, drawn at
%               their MEAN OBSERVED depth, plus the part of each column that no
%               level covers. This is the vertical resolution available to the
%               mode fit -- what build_modes and fit_uv_modes really see.
%
% The two are not the same list and not the same count. The gridded product
% carries RBR Solos and Chi-PODs only: every SBE and every ADCP was used to
% BUILD it (SBE pressure fixes the depth interpolation, SBE conductivity is the
% source of S, the ADCPs are the source of u and v) but none of them is a level
% you can fit a mode to. F6 is the extreme case -- 31 instruments on the wire,
% 8 levels in the product, because its 19 SBE 56s are not in there.
%
% In 'gridded' mode the UNCOVERED part of each column is drawn as a pale band:
% above the shallowest level, and below the deepest. Those are the depths where
% a fitted mode is extrapolating, not interpolating. At F1 that band runs from
% 1380 m to 3390 m -- 59% of the water column -- which is the single most
% important thing this figure has to say.
%
% OPTIONS (all optional; pass a struct, missing fields take the default)
%
%   mode        'deployed' | 'gridded'                       ['deployed']
%   box         [lon0 lon1 lat0 lat1]                        [134 135 8 9] %adjusted
%   zmax        bottom of the box [m]; [] fits the deepest mooring    []
%   bathy       netCDF with lat, lon, z. [] takes cfg.bathy from
%               fleat_config_base, which is where every other plot gets it
%   nlev        number of DISCRETE colours / contour intervals     [12]
%   bathystyle  'contour3' | 'surf' | 'both'                   ['both']
%   view        [az el] starting view                           [-37.5 22]
%   vexag       vertical stretch of the box, purely cosmetic     [0.75]
%   showdead    draw instruments the deployment list flagged as not
%               recording, as hollow grey markers                 [true]
%   moorings    cellstr subset, e.g. {'F2','F4'}            [all five]
%   labels      annotate each string with its name and depth    [true]
%
% ADCP ORIENTATION IS INFERRED, NOT RECORDED. Nothing in the FLEAT dataset
% says which way an ADCP looked. fleat_layout_data applies the rule "a live
% ADCP within 60 m of the recorded bottom faces up, everything else faces
% down", which is consistent with where velocity exists at all five moorings
% and is strongly supported at F4: the deep failure on 2016-11-17 removes
% velocity from 636-823 m, exactly the range an up-looker at 857 m would
% cover. It is still an inference, and this figure draws it as one -- the
% direction carets are open, not filled. Check the mooring diagrams before
% putting an orientation in a paper.
%
% RETURNS a struct of handles: h.fig h.ax h.bathy h.cb h.stem h.inst h.leg,
% plus h.data, the table this was drawn from, so the numbers behind any marker
% can be read off without re-deriving them.
%
% SEE ALSO fleat_layout_data, fleat_config_base, flux_movie
%
% author TEJ, with help from Claude

%% ---- options ---------------------------------------------------------
if nargin < 1 || isempty(opt), opt = struct(); end
if ~isstruct(opt)
    error('pmooring_layout:opt', ...
        ['Takes a single options STRUCT, e.g.\n' ...
         '    pmooring_layout(struct(''mode'',''gridded''))\n' ...
         'Got a %s. The options are listed in the help text.'], class(opt));
end

def = struct('mode','deployed', 'box',[134.3 134.75 8.4 8.8], 'zmax',[], 'bathy',[], ...
             'nlev',12, 'bathystyle','both', 'view',[-37.5 22], 'vexag',0.75, ...
             'showdead',true, 'moorings',[], 'labels',true);
fn = fieldnames(def);
unknown = setdiff(fieldnames(opt), fn);
if ~isempty(unknown)
    % Silently ignoring a misspelled option is how you spend ten minutes
    % wondering why opt.moorings did nothing.
    error('pmooring_layout:unknownopt', ...
        'Unknown option(s): %s.\nKnown options: %s.', ...
        strjoin(unknown', ', '), strjoin(fn', ', '));
end
for k = 1:numel(fn)
    if ~isfield(opt, fn{k}) || isempty(opt.(fn{k}))
        opt.(fn{k}) = def.(fn{k});
    end
end

opt.mode = validatestring(opt.mode, {'deployed','gridded'}, ...
                          'pmooring_layout', 'opt.mode');
opt.bathystyle = validatestring(opt.bathystyle, {'contour3','surf','both'}, ...
                          'pmooring_layout', 'opt.bathystyle');

if numel(opt.box) ~= 4 || opt.box(2) <= opt.box(1) || opt.box(4) <= opt.box(3)
    error('pmooring_layout:box', ...
        ['opt.box is [lon0 lon1 lat0 lat1] with lon0<lon1 and lat0<lat1. ' ...
         'Got %s.'], mat2str(opt.box));
end

D = fleat_layout_data();
all_m = fieldnames(D)';
for k = 1:numel(all_m)
    % Normalise once, so nothing downstream has to care whether cell2table
    % gave back a logical column or a cell of logicals.
    if iscell(D.(all_m{k}).dep.dead)
        D.(all_m{k}).dep.dead = logical(cell2mat(D.(all_m{k}).dep.dead));
    end
end
if isequal(opt.moorings, def.moorings)
    moor = all_m;
else
    moor = cellstr(opt.moorings);
    bad = setdiff(moor, all_m);
    if ~isempty(bad)
        % F5 is the one people reach for, and it does not exist.
        error('pmooring_layout:mooring', ...
            ['Unknown mooring(s): %s. The array is %s -- there is no F5.'], ...
            strjoin(bad, ', '), strjoin(all_m, ' '));
    end
end

%% ---- instrument styling ---------------------------------------------
% One place that decides what every instrument type looks like, so the legend
% and the markers cannot drift apart: both are built from this struct.
%
% The two thermistor types are the loud ones because they are the measurement
% the mode fit is built on. Everything else is deliberately quieter.
S = struct();
S.RBRSolo = style('RBR Solo (T)',       'o', [0.85 0.10 0.10], 5.5);
S.SBE56   = style('SBE 56 (T)',         'o', [1.00 0.84 0.00], 5.5);
S.SBE39   = style('SBE 39 (T, P)',      '^', [0.95 0.55 0.10], 5.5);
S.SBE37   = style('SBE 37 (T, S, P)',   'd', [0.20 0.70 0.90], 6.5);
S.RDIadcp = style('RDI ADCP',           's', [1.00 1.00 1.00], 9.0);
S.RDIadcp.edge  = [0.05 0.25 0.85];     % blue box around a white fill
S.RDIadcp.width = 1.6;
S.MP      = style('McLane profiler',    's', [1.00 0.84 0.00], 6.0);

ADCPBLUE = S.RDIadcp.edge;
DEADGREY = [0.60 0.60 0.60];

%% ---- bathymetry -----------------------------------------------------
if isempty(opt.bathy)
    base = fleat_config_base('F2');     % any mooring; cfg.bathy is array-wide
    opt.bathy = base.bathy;
end
if exist(opt.bathy, 'file') ~= 2
    error('pmooring_layout:nobathy', ...
        ['No bathymetry at\n    %s\nPass one as opt.bathy, or fix cfg.bathy ' ...
         'in fleat_config_base.'], opt.bathy);
end

bl = double(ncread(opt.bathy,'lat'));
bo = double(ncread(opt.bathy,'lon'));
bz = double(ncread(opt.bathy,'z'));         % [lon x lat], positive DOWN

il = bl >= opt.box(3) & bl <= opt.box(4);
io = bo >= opt.box(1) & bo <= opt.box(2);
if ~any(il) || ~any(io)
    error('pmooring_layout:boxoutside', ...
        ['opt.box is outside the bathymetry, which spans lon %.2f-%.2f, ' ...
         'lat %.2f-%.2f.'], min(bo), max(bo), min(bl), max(bl));
end
bl = bl(il);  bo = bo(io);  bz = bz(io, il);

% Decimate. The full grid is ~2000x2000; contour3 on 4M points crawls, and at
% the scale this figure is read the detail is invisible anyway.
dec = max(1, round(numel(bo)/300));
bo = bo(1:dec:end);  bl = bl(1:dec:end);  bz = bz(1:dec:end, 1:dec:end);
BZ = bz.';                                   % [lat x lon], as surf/contour want
fprintf('bathymetry: %d x %d after subsetting, decimated by %d\n', ...
        numel(bo), numel(bl), dec);

%% ---- depth range of the box -----------------------------------------
Hs = cellfun(@(m) D.(m).H, moor);
if isempty(opt.zmax)
    % Round up past the deepest mooring AND past the deepest bathymetry in the
    % box, so nothing is drawn outside the axes.
    deepest = max(BZ(isfinite(BZ)));
    if isempty(deepest), deepest = 0; end
    opt.zmax = 250 * ceil(max([Hs(:); deepest]) / 250);
end
if opt.zmax <= 0
    error('pmooring_layout:zmax', 'opt.zmax is a positive depth [m]. Got %g.', opt.zmax);
end

%% ---- figure ---------------------------------------------------------
h.fig = figure('Color','w', 'Position',[60 60 980 820], ...
               'Name', sprintf('FLEAT array -- %s', opt.mode));
h.ax  = axes(h.fig); hold(h.ax,'on'); box(h.ax,'on'); grid(h.ax,'on');

clev = linspace(0, opt.zmax, opt.nlev+1);
colormap(h.ax, cmocean('deep', opt.nlev));   % discrete: one colour per interval
caxis(h.ax, [0 opt.zmax]);

h.bathy = gobjects(0);
if any(strcmp(opt.bathystyle, {'surf','both'}))
    % The seafloor itself. Flat shading keeps the discrete colormap discrete --
    % 'interp' would blend across the levels and defeat the point of nlev.
    h.bathy(end+1) = surf(h.ax, bo, bl, BZ, 'EdgeColor','none', ...
                          'FaceColor','flat', 'FaceAlpha',0.92);
end
if any(strcmp(opt.bathystyle, {'contour3','both'}))
    % contour3 puts each contour at its OWN depth, so these sit on the
    % seafloor rather than being projected onto the floor of the box. That is
    % what makes the slope readable when the figure is rotated.
    [~, h.bathy(end+1)] = contour3(h.ax, bo, bl, BZ, clev, ...
                                   'LineColor',[0.15 0.15 0.15], 'LineWidth',0.4);
end

h.cb = colorbar(h.ax);
ylabel(h.cb, 'bottom depth (m)');
set(h.cb, 'Direction','reverse', 'Ticks', clev(1:2:end));

%% ---- moorings -------------------------------------------------------
h.stem = gobjects(0);
h.inst = gobjects(0);
rows   = {};                                 % accumulates h.data

for k = 1:numel(moor)
    m  = moor{k};
    Mk = D.(m);
    x  = Mk.lon;  y = Mk.lat;  H = Mk.H;

    % The wire, drawn full depth to the recorded bottom whether or not
    % anything is hanging on the lower part of it.
    h.stem(end+1) = plot3(h.ax, [x x], [y y], [0 H], '-', ...
                          'Color',[0.35 0.35 0.35], 'LineWidth',0.9);
    plot3(h.ax, x, y, H, 'v', 'MarkerFaceColor',[0.35 0.35 0.35], ...
          'MarkerEdgeColor','none', 'MarkerSize',5);

    switch opt.mode
    % ------------------------------------------------------------------
    case 'deployed'
        T = Mk.dep;
        for j = 1:height(T)
            typ = T.type{j};  z = T.z(j);  dead = T.dead(j);
            if dead && ~opt.showdead, continue; end

            if strcmp(typ,'MP')
                % The McLane profiler has no fixed depth -- it crawls the wire.
                % At F1 it was meant to cover the gap between the deepest wire
                % thermistor and the anchor instruments, so it is drawn as a
                % tall box spanning that gap rather than as a point.
                % THE PROFILER'S PROGRAMMED RANGE IS NOT RECORDED ANYWHERE, so
                % this box is not a claim about where it travelled. It spans
                % the biggest gap in the GRIDDED levels, which is the hole the
                % profiler was on the wire to fill and the hole that is still
                % there because it failed. At F1 that is 1380 -> 3362 m.
                %
                % Deliberately NOT the gap in the deployed list: that would
                % come out as 2304 -> 3362 m, because a live SBE 37 sits at
                % 2304 m. It recorded, but it is not a level in the gridded
                % product, so it does nothing to close the gap this box is
                % about.
                zf = sort(Mk.grd.z(isfinite(Mk.grd.z)));
                if numel(zf) < 2, continue; end
                [gap, ig] = max(diff(zf));
                if gap < 0.05*Mk.H, continue; end   % no gap worth drawing
                zt = zf(ig);  zb = zf(ig+1);
                h.inst(end+1) = mpbox(h.ax, x, y, zt, zb, opt, S.MP.face, dead);
                rows(end+1,:) = {m, typ, T.sn{j}, zt, zb, dead, T.note{j}, ''};
                continue
            end

            if ~isfinite(z), continue; end
            if ~isfield(S, typ)
                error('pmooring_layout:type', ...
                    ['No style for instrument type "%s" (%s, SN %s). Add one ' ...
                     'to the S struct in this file -- the legend is built from ' ...
                     'it, so a type with no style would also be missing from ' ...
                     'the key.'], typ, m, T.sn{j});
            end
            st = S.(typ);
            if dead
                % Hollow and grey: present on the wire, absent from the record.
                h.inst(end+1) = plot3(h.ax, x, y, z, st.mark, ...
                    'MarkerFaceColor','none', 'MarkerEdgeColor',DEADGREY, ...
                    'MarkerSize',st.size, 'LineWidth',0.8);
            else
                h.inst(end+1) = plot3(h.ax, x, y, z, st.mark, ...
                    'MarkerFaceColor',st.face, 'MarkerEdgeColor',st.edge, ...
                    'MarkerSize',st.size, 'LineWidth',st.width);
            end

            % ADCP look direction, as an OPEN caret offset just off the wire.
            % Open because the orientation is inferred -- see the header.
            if strcmp(typ,'RDIadcp') && ~dead && ~isempty(T.orient{j})
                dx = 0.016 * (opt.box(2)-opt.box(1));
                if strcmp(T.orient{j},'up'), mk = '^'; dz = -0.018*opt.zmax;
                else,                        mk = 'v'; dz =  0.018*opt.zmax;
                end
                plot3(h.ax, x+dx, y, z+dz, mk, 'MarkerFaceColor','none', ...
                      'MarkerEdgeColor',ADCPBLUE, 'MarkerSize',7, ...
                      'LineWidth',1.2);
                plot3(h.ax, [x x+dx], [y y], [z z+dz], ':', ...
                      'Color',ADCPBLUE, 'LineWidth',0.8);
            end

            rows(end+1,:) = {m, typ, T.sn{j}, z, z, dead, T.note{j}, T.orient{j}};
        end

    % ------------------------------------------------------------------
    case 'gridded'
        T  = Mk.grd;
        zz = T.z;
        % THE FITTABLE SPAN IS SET BY VELOCITY, not by where a thermistor
        % happens to hang. F1 has two thermistors on the anchor at 3362 and
        % 3374 m carrying temperature only; counting them would make the
        % column look covered to 3374 m and hide the real hole, which runs
        % from the deepest velocity level at 1380 m to the 3390 m bottom --
        % 59% of the water column, and the most important thing here.
        zu = zz(T.covU > 0);
        if isempty(zu), zu = zz; end
        % A LEVEL IS DEFICIENT RELATIVE TO ITS OWN MOORING, not to 100%. An
        % absolute threshold flags nearly every level at F1 (whose whole record
        % stops on 2016-12-10, so every level reads ~60%) and at F6 (deployed
        % 18 days late, so every level reads ~94%) -- which says something true
        % about the mooring and nothing about the level. Comparing each level
        % with the mooring's own median separates the two.
        ref = median(T.covU(T.covU > 0));
        if isempty(ref) || ~isfinite(ref), ref = 100; end
        thin = T.covU > 0 & T.covU < 0.9*ref;
        band(h.ax, x, y, opt, 0,       min(zu), [0.85 0.33 0.10]);
        band(h.ax, x, y, opt, max(zu), H,       [0.85 0.33 0.10]);

        for j = 1:height(T)
            % A short horizontal bar per level: this is a LEVEL, not a point
            % instrument, and it is what a mode gets fitted to.
            dx = 0.013 * (opt.box(2)-opt.box(1));
            if T.covU(j) == 0
                % Temperature only. Real data, but no velocity ever reached it,
                % so it cannot carry a u/v mode fit. Thin and grey.
                h.inst(end+1) = plot3(h.ax, [x-dx x+dx], [y y], [zz(j) zz(j)], '-', ...
                    'Color',[0.55 0.55 0.55], 'LineWidth',0.9);
            else
                h.inst(end+1) = plot3(h.ax, [x-dx x+dx], [y y], [zz(j) zz(j)], '-', ...
                    'Color',[0.10 0.35 0.75], 'LineWidth',1.6);
            end
            if thin(j)
                % Velocity does not reach this level for the whole record, so
                % it is not a level the flux can lean on. Ring it.
                plot3(h.ax, x, y, zz(j), 'o', 'MarkerEdgeColor',[0.85 0.10 0.10], ...
                      'MarkerFaceColor','none', 'MarkerSize',5, ...
                      'LineWidth',1.0);
            end
            rows(end+1,:) = {m, 'gridlevel', T.sn{j}, zz(j), zz(j), false, ...
                             sprintf('T %.0f%% S %.0f%% u %.0f%%', ...
                                     T.covT(j), T.covS(j), T.covU(j)), ''};
        end
    end

    if opt.labels
        text(h.ax, x, y, -0.035*opt.zmax, sprintf('%s', m), ...
             'FontWeight','bold', 'FontSize',11, 'HorizontalAlignment','center');
        if strcmp(opt.mode,'gridded')
            zu = Mk.grd.z(Mk.grd.covU > 0);
            if isempty(zu), zu = Mk.grd.z; end
            lab = sprintf('%d levels, %d with u\\newlineu: %.0f-%.0f m\\newlineH %d m', ...
                          height(Mk.grd), numel(zu), min(zu), max(zu), H);
        else
            live = ~Mk.dep.dead;
            lab = sprintf('%d instr (%d live)\\newlineH %d m', ...
                          height(Mk.dep), sum(live), H);
        end
        text(h.ax, x, y, 0.045*opt.zmax, lab, 'FontSize',8, ...
             'Color',[0.25 0.25 0.25], 'HorizontalAlignment','center', ...
             'VerticalAlignment','top');
    end
end

h.data = cell2table(rows, 'VariableNames', ...
    {'mooring','type','sn','ztop','zbot','dead','note','orient'});

%% ---- axes -----------------------------------------------------------
xlabel(h.ax, 'longitude (\circE)');
ylabel(h.ax, 'latitude (\circN)');
zlabel(h.ax, 'depth (m)');
xlim(h.ax, opt.box(1:2));
ylim(h.ax, opt.box(3:4));
zlim(h.ax, [-0.06*opt.zmax opt.zmax]);
set(h.ax, 'ZDir','reverse', 'Layer','top', 'FontSize',10);
view(h.ax, opt.view);

% NOT a physical aspect ratio, and it cannot be. One degree is ~111 km against
% a 3.4 km water column, so true scaling would squash every mooring into the
% seafloor. pbaspect sets how the BOX looks; opt.vexag stretches the vertical.
pbaspect(h.ax, [1 1 opt.vexag]);
rotate3d(h.ax, 'on');

if strcmp(opt.mode,'gridded')
    ttl = 'FLEAT array: gridded levels available to the mode fit';
    sub = ['bars are levels in F*\_fullfield.mat at mean observed depth; ' ...
           'red rings mark u coverage well below that mooring''s own ' ...
           'median; pale bands are column with ' ...
           'no velocity level, where a fit extrapolates'];
else
    ttl = 'FLEAT array: instruments as deployed';;
end
title(h.ax, {ttl, ['\rm\fontsize{9}']});

%% ---- legend ---------------------------------------------------------
% Built from off-screen dummies so the legend lists each type once, in a fixed
% order, however many of each happen to be in the chosen subset.
lh = gobjects(0); lt = {};
if strcmp(opt.mode,'deployed')
    order = {'RBRSolo','SBE56','SBE39','SBE37','RDIadcp','MP'};
    present = unique(h.data.type);
    for k = 1:numel(order)
        if ~ismember(order{k}, present), continue; end
        st = S.(order{k});
        lh(end+1) = plot3(h.ax, NaN, NaN, NaN, st.mark, ...
            'MarkerFaceColor',st.face, 'MarkerEdgeColor',st.edge, ...
            'MarkerSize',st.size, 'LineWidth',st.width);
        nm = st.name;
        if strcmp(order{k},'MP') && all(h.data.dead(strcmp(h.data.type,'MP')))
            nm = [nm ' (no data)'];
        end
        lt{end+1} = nm;
    end
    if opt.showdead && any(h.data.dead)
        lh(end+1) = plot3(h.ax, NaN,NaN,NaN, 'o', 'MarkerFaceColor','none', ...
            'MarkerEdgeColor',DEADGREY, 'MarkerSize',5.5, 'LineWidth',0.8);
        lt{end+1} = 'did not record';
    end
    lh(end+1) = plot3(h.ax, NaN,NaN,NaN, '^', 'MarkerFaceColor','none', ...
        'MarkerEdgeColor',ADCPBLUE, 'MarkerSize',7, 'LineWidth',1.2);
    lt{end+1} = 'ADCP looks up / down (inferred)';
else
    lh(end+1) = plot3(h.ax, NaN,NaN,NaN, '-', 'Color',[0.10 0.35 0.75], 'LineWidth',1.6);
    lt{end+1} = 'level with u (mode fit basis)';
    lh(end+1) = plot3(h.ax, NaN,NaN,NaN, '-', 'Color',[0.55 0.55 0.55], 'LineWidth',0.9);
    lt{end+1} = 'level with T only (no u)';
    lh(end+1) = plot3(h.ax, NaN,NaN,NaN, 'o', 'MarkerEdgeColor',[0.85 0.10 0.10], ...
        'MarkerFaceColor','none', 'MarkerSize',5, 'LineWidth',1.0);
    lt{end+1} = 'u coverage < 90% of mooring median';
    lh(end+1) = patch(h.ax, 'XData',[NaN NaN NaN], 'YData',[NaN NaN NaN], ...
        'ZData',[NaN NaN NaN], 'FaceColor',[0.85 0.33 0.10], 'FaceAlpha',0.18, ...
        'EdgeColor','none');
    lt{end+1} = 'no level: fit extrapolates';
end
h.leg = legend(h.ax, lh, lt, 'Location','northeastoutside', 'FontSize',9, ...
               'AutoUpdate','off');

%% ---- what got drawn -------------------------------------------------
fprintf('\n%s mode, %d moorings, box lon %.2f-%.2f lat %.2f-%.2f, z to %g m\n', ...
        opt.mode, numel(moor), opt.box, opt.zmax);
for k = 1:numel(moor)
    m = moor{k};
    if strcmp(opt.mode,'gridded')
        g  = D.(m).grd;
        zu = g.z(g.covU > 0);
        if isempty(zu), zu = g.z; end
        ref = median(g.covU(g.covU > 0));
        if isempty(ref) || ~isfinite(ref), ref = 100; end
        fprintf(['  %s: %2d levels, u over %6.1f-%6.1f m of H=%4d m ' ...
                 '(%4.1f%% of column), %d thin in u\n'], ...
                m, height(g), min(zu), max(zu), D.(m).H, ...
                100*(max(zu)-min(zu))/D.(m).H, sum(g.covU > 0 & g.covU < 0.9*ref));
    else
        t = D.(m).dep;
        fprintf('  %s: %2d instruments, %2d live, %2d dead, %d ADCP (%d up, %d down)\n', ...
                m, height(t), sum(~t.dead), sum(t.dead), ...
                sum(strcmp(t.type,'RDIadcp')), sum(strcmp(t.orient,'up')), ...
                sum(strcmp(t.orient,'down')));
    end
end
fprintf('rotate3d is on -- drag to spin. h.data holds the numbers behind the markers.\n\n');

end


%% ======================================================================
function s = style(name, mark, face, sz)
% One instrument type's appearance. Defaults that most types share live here
% so each line above only states what is different about it.
s = struct('name',name, 'mark',mark, 'face',face, 'size',sz, ...
           'edge',[0.15 0.15 0.15], 'width',0.5);
end


function ph = mpbox(ax, x, y, zt, zb, opt, col, dead)
% MPBOX  The McLane profiler as a tall box spanning the depths it worked.
% A point marker would be a lie: the instrument's whole purpose is that it
% moves, and at F1 it is the only thing that was ever going to sample the
% 1380-3360 m gap.
dx = 0.017 * (opt.box(2)-opt.box(1));
dy = 0.017 * (opt.box(4)-opt.box(3));
X = [x-dx x+dx x+dx x-dx];
Y = [y-dy y-dy y+dy y+dy];
% BOTH F1 PROFILERS RETURNED NO DATA, so in practice this box is always the
% dead case. It keeps the yellow anyway -- the box is there to show the depth
% range the array was MEANT to cover, and greying it out would hide the gap
% instead of pointing at it. The dashed edge and the pale fill are what say no
% data came back.
if dead
    fc = col; al = 0.22; ec = [0.55 0.45 0.00]; ls = '--';
else
    fc = col; al = 0.55; ec = [0.35 0.30 0.00]; ls = '-';
end
% Two facing panels rather than a closed box: a solid box hides the wire and
% the instruments behind it from most viewing angles.
%
% At F1 both profilers (MM01, MM02) returned no data, which is why the
% deployment list has them commented out. The box still spans the gap, because
% the point of drawing it is that the gap was meant to be covered and is not.
ph = patch(ax, 'XData',X, 'YData',Y, 'ZData',[zt zt zt zt], ...
           'FaceColor',fc, 'FaceAlpha',al, 'EdgeColor',ec, 'LineStyle',ls);
patch(ax, 'XData',[x-dx x+dx x+dx x-dx], 'YData',[y y y y], ...
      'ZData',[zt zt zb zb], 'FaceColor',fc, 'FaceAlpha',al, ...
      'EdgeColor',ec, 'LineStyle',ls);
patch(ax, 'XData',[x x x x], 'YData',[y-dy y+dy y+dy y-dy], ...
      'ZData',[zt zt zb zb], 'FaceColor',fc, 'FaceAlpha',al, ...
      'EdgeColor',ec, 'LineStyle',ls);
end


function band(ax, x, y, opt, z1, z2, col)
% BAND  A pale panel over a depth range no gridded level covers, i.e. where a
% fitted mode is extrapolating. Drawn as two crossed panels so it stays
% visible from any rotation.
if ~isfinite(z1) || ~isfinite(z2) || z2 - z1 < 1e-6, return; end
dx = 0.020 * (opt.box(2)-opt.box(1));
dy = 0.020 * (opt.box(4)-opt.box(3));
patch(ax, 'XData',[x-dx x+dx x+dx x-dx], 'YData',[y y y y], ...
      'ZData',[z1 z1 z2 z2], 'FaceColor',col, 'FaceAlpha',0.18, 'EdgeColor','none');
patch(ax, 'XData',[x x x x], 'YData',[y-dy y+dy y+dy y-dy], ...
      'ZData',[z1 z1 z2 z2], 'FaceColor',col, 'FaceAlpha',0.18, 'EdgeColor','none');
end
