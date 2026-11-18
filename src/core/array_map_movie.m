function out = array_map_movie(R, mopt)
% ARRAY_MAP_MOVIE  The FLEAT array on a bathymetry map: one vector per mooring,
%                  optionally one scalar dot per mooring, as a frame or a movie.
%
%   out = array_map_movie(R)                     % static, record-mean flux
%   out = array_map_movie(R, mopt)
%
%   R     cell array (or struct array) of FLEAT_IT_FLUX outputs, exactly as
%         built by section 1 of turner_fleat_IT_flux. Failed moorings may be in
%         it; they are carried as a marker with no vector, never dropped
%         silently.
%
%   out .A        the aligned array struct (ALIGN_FLUX_RECORDS)
%       .V        [nm x nt]  what was drawn: Vx, Vy, mag, units, label
%       .S        the scalar series and its colour limits, [] if none asked for
%       .dot      [nm x 2]   where each scalar dot was placed, lon/lat
%       .frames   the frame indices animated
%       .outfile  the mp4 written, '' if none
%
%   ---------------------------------------------------------------------
%   WHAT IT DRAWS
%   ---------------------------------------------------------------------
%   mopt.vector picks the arrow. All five are per-mooring time series, so all
%   five animate; the units and the reference arrow change with the choice.
%
%     'flux'       (default) depth-INTEGRATED baroclinic M2 energy flux, kW/m.
%                  This is out.flux.Fx/Fy, the quantity Alford and Nash map and
%                  the one the pipeline exists to produce.
%     'flux_mean'  the same divided by cfg.H, W/m^2, i.e. the depth AVERAGE.
%                  Offered because "depth-average flux" is the natural phrase
%                  and is NOT what 'flux' is. Across this array the distinction
%                  is not cosmetic: H runs 434 m at F6 to 3390 m at F1, so
%                  dividing by it reorders the moorings. Say which one a figure
%                  shows.
%     'bt'         barotropic M2 velocity, mode 0 of the velocity fit, m/s.
%     'uv_slice'   baroclinic u',v' at mopt.z metres, m/s, off U.up/U.vp.
%     'uv_mean'    depth mean of the OBSERVED hourly u,v over the sampled
%                  range, m/s. Read the note on this one below -- it is not the
%                  depth mean of u'.
%
%   mopt.scalar adds a filled dot beside each mooring, coloured by T, S, chi or
%   eps at that frame. It has its own colour scale and its own colourbar, and
%   the dots are offset so that they collide neither with each other nor with
%   the arrows (see LAYOUT below).
%
%   ---------------------------------------------------------------------
%   WHY 'uv_mean' IS THE OBSERVED VELOCITY AND NOT u'
%   ---------------------------------------------------------------------
%   The depth mean of the baroclinic velocity is ZERO. Not approximately: the
%   baroclinic modes Z^u are defined with zero depth integral (Nash eq. 6, and
%   fit_uv_modes checks it -- U.bt_check is 6.7e-5 m/s at F6), so
%   trapz(z, U.up)/H is zero by construction and a map of it would be a map of
%   numerical noise. Anyone asking for "depth-mean velocity" wants the real
%   flow, so 'uv_mean' averages the observed Mh.u and Mh.v over the depths that
%   actually reported in that hour.
%
%   That makes it a different kind of quantity from the other four: it is the
%   total flow, mean plus mesoscale plus barotropic tide plus whatever the
%   baroclinic field leaves after an incomplete-depth average, and it is NOT
%   bandpassed. Use 'bt' for the barotropic tide alone. The average is also
%   over the SAMPLED range, which at F6 is 137-420 m of 434 m, so it is not the
%   depth mean of the column either. The axis label says "sampled range" for
%   that reason.
%
%   ---------------------------------------------------------------------
%   BATHYMETRY: DISCRETE BANDS, NON-UNIFORM LEVELS
%   ---------------------------------------------------------------------
%   FLUX_MOVIE contours 0:250:6000, which is uniform and wrong for this box.
%   The array sits on the Velasco Reef flank: within the requested extent the
%   grid runs from reef-flat depths to beyond 2000 m, and F6's own note records
%   15 m to 1578 m inside two kilometres. A uniform 250 m contour puts every
%   mooring but F1 in one or two bands and renders the flank as a single cliff.
%
%   So the levels are non-uniform and set by where the moorings are: tight over
%   the reef and the shelf break, loose in the basin. Getting MATLAB to colour
%   non-uniform levels DISCRETELY takes one trick, because contourf assigns
%   colour by VALUE, so unequal bands get unequal colour steps and the map
%   reads as a smooth ramp again. Instead we contour the BAND INDEX --
%   discretize(z, levels) -- on unit levels, which gives every band exactly one
%   colour of equal visual weight, and then relabel the colourbar ticks with
%   the depths those bands mean. The colourbar is therefore non-linear in depth
%   on purpose, and it is tick-labelled with the real numbers so that is
%   visible rather than hidden.
%
%   ---------------------------------------------------------------------
%   LAYOUT: WHY THE DOTS DO NOT MOVE
%   ---------------------------------------------------------------------
%   F4 and F6 are 1.7 km apart -- 0.015 deg of latitude, under 4% of the height
%   of this box -- and F2, F3, F4, F6 all sit inside 0.06 deg. A scalar dot
%   drawn at the mooring position would overlap its neighbour's and sit under
%   its own arrow.
%
%   Each dot is placed once, on a small circle around its mooring, at the angle
%   that maximises its distance from every other mooring, every dot already
%   placed, and the local arrow. Crowded moorings choose first, since they have
%   the least room. Distances are measured in aspect-corrected degrees
%   (lon*cosd(lat)), so "distance" on screen is distance in the metric.
%
%   The layout is computed ONCE, from the record-mean arrow directions, and
%   held fixed for every frame. Recomputing per frame would be more optimal per
%   frame and unreadable as a movie: the dots would jump around the markers as
%   the arrows rotated, and a viewer would read the motion as data. A fixed dot
%   that an arrow occasionally crosses is the better trade.
%
%   ---------------------------------------------------------------------
%   WHAT A FRAME HAS TO ADMIT
%   ---------------------------------------------------------------------
%   One scale for all moorings, fixed for the whole record, set from the 99th
%   percentile of the pooled magnitude, with the reference arrow and its units
%   in the corner of every frame. Otherwise frames are not comparable to each
%   other and moorings are not comparable within a frame.
%
%   A mooring with no data in a frame gets its marker and a greyed label and NO
%   arrow -- not a zero-length one. "No data" and "no flux" are different
%   claims and the ALIGN_FLUX_RECORDS NaN padding exists to keep them apart.
%
%
% author TEJ, with help from Claude
% 30 October 2026

%% ---- options ---------------------------------------------------------
if nargin < 2 || isempty(mopt), mopt = struct(); end
mopt = setdef(mopt, 'vector',    'flux');
mopt = setdef(mopt, 'scalar',    '');
mopt = setdef(mopt, 'z',         []);        % depth for 'uv_slice'
mopt = setdef(mopt, 'scalar_z',  []);        % depth for the scalar; [] = mean
mopt = setdef(mopt, 'scalar_log', []);       % [] = decide from the variable
mopt = setdef(mopt, 'box',       [134.3 134.75 8.4 8.8]);
mopt = setdef(mopt, 'levels',    []);        % [] = the default band set
mopt = setdef(mopt, 'bathy',     '');
mopt = setdef(mopt, 'frames',    []);        % [] = one static record-mean frame
mopt = setdef(mopt, 'stride_hr', 6);         % frame spacing when animating
mopt = setdef(mopt, 'scale',     []);        % deg per unit; [] = auto
mopt = setdef(mopt, 'dot_r',     0.030);     % dot offset radius, deg of lat
mopt = setdef(mopt, 'dot_size',  150);
mopt = setdef(mopt, 'isobath',   false);     % draw each mooring's own isobath
mopt = setdef(mopt, 'outfile',   '');        % '' = animate without writing
mopt = setdef(mopt, 'fps',       10);
mopt = setdef(mopt, 'series',    true);      % the time-series panel with cursor
mopt = setdef(mopt, 'dt_hr',     1);

if isstruct(R), R = num2cell(R); end
if ~iscell(R) || isempty(R)
    error('array_map_movie:input', ...
        'Pass the cell array of fleat_IT_flux outputs, e.g. array_map_movie(R).');
end

%% ---- the common clock and the positions ------------------------------
% ALIGN_FLUX_RECORDS already owns this: the union grid, the NaN padding, the
% snap-not-interpolate rule and the per-mooring coverage. Everything else this
% function needs is put on A.time by SNAP below, so there is one clock.
A  = align_flux_records(R, mopt.dt_hr);
nm = numel(A.name);
nt = numel(A.time);

%% ---- the vector field ------------------------------------------------
V = get_vector(R, A, mopt);
fprintf('vector: %s [%s], %d of %d moorings carry it\n', ...
        V.label, V.units, nnz(any(isfinite(V.Vx),2)), nm);

%% ---- the scalar field, if asked for ----------------------------------
S = [];
if ~isempty(mopt.scalar)
    S = get_scalar(R, A, mopt);
    fprintf('scalar: %s, colour limits %.3g to %.3g%s\n', ...
            S.label, S.clim(1), S.clim(2), tern(S.log, ' (log10)', ''));
end

%% ---- bathymetry ------------------------------------------------------
B = get_bathy(R, mopt);

%% ---- arrow scale -----------------------------------------------------
% One number for the whole record and the whole array. The 99th percentile, not
% the max, so a single outlying hour cannot shrink every other frame to nothing
% -- and it is reported, because a reader of the figure cannot recover it.
if isempty(mopt.scale) || ~isfinite(mopt.scale)
    ref = prctile(V.mag(:), 99);
    if ~isfinite(ref) || ref <= 0, ref = 1; end
    V.ref   = round_nice(ref);
    mopt.scale = 0.18*diff(mopt.box(3:4)) / V.ref;
else
    V.ref = round_nice(prctile(V.mag(:), 99));
end
fprintf('arrow scale: %.4g deg per %s (99th pct %.3g, reference arrow %.3g)\n', ...
        mopt.scale, V.units, prctile(V.mag(:),99), V.ref);

%% ---- frames ----------------------------------------------------------
% [] means one static frame showing the record mean, which is the figure that
% should exist before any movie does. Anything else animates.
static = isempty(mopt.frames);
if static
    frames = NaN;                     % sentinel: the mean, not an index
else
    if ischar(mopt.frames) || isstring(mopt.frames)
        step   = max(1, round(mopt.stride_hr/mopt.dt_hr));
        anyfin = any(isfinite(V.Vx), 1);
        k      = find(anyfin);
        if isempty(k)
            error('array_map_movie:noframes', ...
                'No timestep has a finite %s at any mooring.', V.label);
        end
        frames = k(1):step:k(end);
    else
        frames = mopt.frames(:).';
    end
    frames = frames(frames >= 1 & frames <= nt);
    if isempty(frames)
        error('array_map_movie:noframes', ...
            'No requested frame falls inside the %d-step common clock.', nt);
    end
    if numel(frames) > 1, fstep = median(diff(frames)); else, fstep = 1; end
    fprintf('%d frames, every %.0f h, %.1f s at %d fps\n', ...
            numel(frames), mopt.dt_hr*fstep, numel(frames)/mopt.fps, mopt.fps);
end

%% ---- dot layout ------------------------------------------------------
% Fixed for the whole movie. See LAYOUT in the header.
dot = nan(nm,2);
if ~isempty(S)
    mdir = atan2(mean(V.Vy,2,'omitnan'), mean(V.Vx,2,'omitnan'));
    dot  = place_dots(A.lon, A.lat, mdir, mopt);
end

%% ---- the figure ------------------------------------------------------
lat0 = mean(mopt.box(3:4));
fh = figure('Color','w','Position',[60 60 820 tern(mopt.series, 860, 700)]);
% Locked geometry: a figure resized mid-render changes the captured frame
% size, which shears the movie. See capture_frame.
set(fh, 'Resize', 'off', 'PaperPositionMode', 'auto');
FRAME_DPI = 150;   % capture resolution, dots per inch
fsz = [];          % frame size, locked by the first capture_frame call
if mopt.series
    ax  = axes(fh,'Position',[0.10 0.38 0.74 0.57]);
    axs = axes(fh,'Position',[0.10 0.08 0.74 0.20]);
else
    ax  = axes(fh,'Position',[0.10 0.11 0.74 0.82]);
    axs = [];
end
hold(ax,'on'); box(ax,'on');

% ---- discrete bathymetry, non-uniform levels ----
nb = numel(B.lev) - 1;
contourf(ax, B.lon, B.lat, B.idx, 0.5:1:(nb+0.5), 'LineStyle','none');
colormap(ax, pick_cmap('deep', nb));
clim(ax, [0.5 nb+0.5]);
ap = get(ax, 'Position');          % both colourbars key off the map's box
cb = colorbar(ax, 'Position', [ap(1)+ap(3)+0.025 ap(2) 0.025 ap(4)]);
cb.Ticks      = 0.5:1:(nb+0.5);
cb.TickLabels = compose('%g', B.lev);
cb.Direction  = 'reverse';
ylabel(cb, 'depth [m]  (bands are not equal)');

axis(ax, mopt.box);
daspect(ax, [1 cosd(lat0) 1]);
xlabel(ax, 'longitude [\circE]'); ylabel(ax, 'latitude [\circN]');

% Each mooring's own isobath, so an arrow can be read against the depth it sits
% on. Off by default: five of them in a box this size is clutter.
if mopt.isobath
    for i = 1:nm
        if isfinite(A.H(i))
            contour(ax, B.lon, B.lat, B.z, A.H(i)*[1 1], 'w-', 'LineWidth', 0.6);
        end
    end
end

% ---- scalar dots, with their own colour scale ----
% Drawn with explicit per-point RGB rather than through the axes colormap,
% because an axes has ONE colormap and the bathymetry already owns it. The
% colourbar beside them is therefore hand-built from the same map, which keeps
% the dots and their legend provably the same mapping.
hs = [];
if ~isempty(S)
    % Leader lines FIRST, so the dot covers the end of its own line instead of
    % the line being drawn across the dot face.
    for i = 1:nm
        if isfinite(dot(i,1))
            plot(ax, [A.lon(i) dot(i,1)], [A.lat(i) dot(i,2)], '-', ...
                 'Color',[0.45 0.45 0.45], 'LineWidth', 0.5);
        end
    end
    hs = scatter(ax, dot(:,1), dot(:,2), mopt.dot_size, ...
                 repmat([.85 .85 .85], nm, 1), 'filled', ...
                 'MarkerEdgeColor','k', 'LineWidth', 0.6);

    % The strip under the map, positioned off the map's own box rather than
    % written as literals: with mopt.series = false the map fills the figure
    % and a hardcoded y would land the colourbar in the middle of it.
    axc = axes(fh, 'Position', [ap(1) ap(2)-0.055 ap(3) 0.016]);
    imagesc(axc, linspace(S.clim(1), S.clim(2), 256), [0 1], 1:256);
    colormap(axc, S.cmap);
    set(axc, 'YTick', [], 'TickDir','out', 'Layer','top', 'YDir','normal');
    xlabel(axc, S.label, 'FontSize', 9);
    clim(axc, [1 256]);
end

% ---- markers, labels, reference arrow, arrows ----
plot(ax, A.lon, A.lat, 'o', 'MarkerFaceColor','w', 'MarkerEdgeColor','k', ...
     'MarkerSize', 6, 'LineWidth', 1.1);
htx = gobjects(nm,1);
for i = 1:nm
    htx(i) = text(ax, A.lon(i)+0.008, A.lat(i), lbl(R, A, i), ...
                  'FontSize', 9, 'FontWeight','bold', 'Color','k');
end

xr = mopt.box(1) + 0.07*diff(mopt.box(1:2));
yr = mopt.box(3) + 0.07*diff(mopt.box(3:4));
quiver(ax, xr, yr, V.ref*mopt.scale/cosd(lat0), 0, 0, 'k', ...
       'LineWidth', 1.5, 'MaxHeadSize', 2);
text(ax, xr, yr - 0.045*diff(mopt.box(3:4)), ...
     sprintf('%.3g %s', V.ref, V.units), 'FontSize', 9);

hq = quiver(ax, A.lon, A.lat, nan(1,nm), nan(1,nm), 0, ...
            'Color',[0.85 0.1 0.1], 'LineWidth', 2.2, 'MaxHeadSize', 0.9);
ht = title(ax, '', 'FontWeight','normal');

% The footnote travels with the figure, not with the person who made it.
un = unverified_note(R, A);
if ~isempty(un)
    annotation(fh, 'textbox', [0.10 0.955 0.78 0.04], 'String', un, ...
        'EdgeColor','none', 'FontSize', 8, 'Color',[0.6 0 0], ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle');
end

% ---- the time-series panel ----
hcur = [];
if mopt.series
    hold(axs,'on'); box(axs,'on')
    co = lines(max(nm,7));
    for i = 1:nm
        plot(axs, A.time, V.mag(i,:), '-', 'Color', co(i,:), 'LineWidth', 0.7);
    end
    xlim(axs, [A.time(1) A.time(end)]);
    datetick(axs, 'x', 'mmm', 'keeplimits');
    ylabel(axs, sprintf('|%s| [%s]', V.short, V.units), 'FontSize', 9);
    legend(axs, A.name, 'Box','off', 'Orientation','horizontal', ...
           'Location','northoutside', 'FontSize', 8);
    hcur = xline(axs, A.time(1), 'k-', 'LineWidth', 1.2);
    if static && ~isempty(hcur), set(hcur, 'Visible', 'off'); end
end

%% ---- draw ------------------------------------------------------------
vw = [];
if ~static && ~isempty(mopt.outfile)
    vw = VideoWriter(mopt.outfile, 'MPEG-4');  vw.FrameRate = mopt.fps;  open(vw);
end

for k = 1:numel(frames)
    it = frames(k);
    % NOT written with tern(): it is a function, so BOTH branches evaluate and
    % S.val on an empty S is an error. A ternary that always evaluates both
    % sides is only safe on values, never on indexing.
    sv = [];
    if static
        vx = mean(V.Vx,2,'omitnan').';   vy = mean(V.Vy,2,'omitnan').';
        if ~isempty(S), sv = mean(S.val, 2, 'omitnan'); end
        tstr = sprintf('record mean, %s to %s', ...
                       datestr(A.time(1),'dd mmm yyyy'), ...
                       datestr(A.time(end),'dd mmm yyyy'));
    else
        vx = V.Vx(:,it).';   vy = V.Vy(:,it).';
        if ~isempty(S), sv = S.val(:,it); end
        tstr = datestr(A.time(it), 'dd mmm yyyy HH:MM');
        if ~isempty(hcur), set(hcur, 'Value', A.time(it)); end
    end

    % Missing moorings get NaN, which quiver draws as nothing. A zero-length
    % arrow would read as "no flux" instead of "no data".
    set(hq, 'UData', vx*mopt.scale/cosd(lat0), 'VData', vy*mopt.scale);

    if ~isempty(S)
        set(hs, 'CData', val2rgb(sv, S.clim, S.cmap));
    end
    for i = 1:nm
        set(htx(i), 'Color', tern(isfinite(vx(i)), [0 0 0], [0.55 0.55 0.55]));
    end
    set(ht, 'String', sprintf('%s      %s      n = %d of %d reporting', ...
            tstr, V.label, nnz(isfinite(vx)), nm));
    % Plain drawnow, not 'limitrate': limitrate may skip the update, and
    % the capture below would then write the previous frame's state.
    drawnow

    if ~isempty(vw)
        [im, fsz] = capture_frame(fh, FRAME_DPI, fsz);
        writeVideo(vw, im);
    end
end

if ~isempty(vw)
    close(vw);
    fprintf('wrote %s\n', mopt.outfile);
end

%% ---- out -------------------------------------------------------------
out = struct('A', A, 'V', V, 'S', S, 'dot', dot, ...
             'frames', tern(static, [], frames), 'outfile', mopt.outfile, ...
             'bathy', B, 'mopt', mopt, 'fig', fh);
end

%% ======================================================================
function V = get_vector(R, A, mopt)
% The arrow, on A.time, in one consistent unit. Each branch states its own
% units and its own label; nothing downstream guesses.
nm = numel(A.name);  nt = numel(A.time);
V  = struct('Vx', nan(nm,nt), 'Vy', nan(nm,nt));

switch lower(mopt.vector)
  case 'flux'
    V.Vx = A.Fx/1e3;  V.Vy = A.Fy/1e3;
    V.units = 'kW m^{-1}';  V.short = 'F';
    V.label = 'depth-integrated M_2 baroclinic energy flux';

  case 'flux_mean'
    % Divide by the water depth, NOT by the sampled range: the modal integral
    % that produced Fx already spans 0 to H, so H is the consistent divisor.
    H = A.H(:);
    V.Vx = A.Fx ./ H;  V.Vy = A.Fy ./ H;
    V.units = 'W m^{-2}';  V.short = 'F/H';
    V.label = 'depth-AVERAGED M_2 baroclinic energy flux';

  case 'bt'
    for i = 1:nm
        U = field_or(R{i}, 'uv');
        if isempty(U) || ~isfield(U,'u_bt'), continue; end
        V.Vx(i,:) = snap(U.u_bt, R{i}.time, A.time);
        V.Vy(i,:) = snap(U.v_bt, R{i}.time, A.time);
    end
    need_full(V.Vx, 'bt', 'out.uv');
    V.units = 'm s^{-1}';  V.short = 'u_0';
    V.label = 'barotropic M_2 velocity (mode 0)';

  case 'uv_slice'
    if isempty(mopt.z)
        error('array_map_movie:noz', ...
            'mopt.vector = ''uv_slice'' needs mopt.z, the depth in metres.');
    end
    for i = 1:nm
        U = field_or(R{i}, 'uv');  M = field_or(R{i}, 'modes');
        if isempty(U) || isempty(M) || ~isfield(U,'up'), continue; end
        % Outside the mode grid there is nothing to report. Clamping would
        % quietly hand back the bottom value for a 400 m request at F6.
        if mopt.z < M.z(1) || mopt.z > M.z(end)
            warning('array_map_movie:zrange', ...
                '%s: mopt.z = %g m is outside its mode grid (%g to %g m); left NaN.', ...
                A.name{i}, mopt.z, M.z(1), M.z(end));
            continue
        end
        V.Vx(i,:) = snap(interp1(M.z, U.up, mopt.z), R{i}.time, A.time);
        V.Vy(i,:) = snap(interp1(M.z, U.vp, mopt.z), R{i}.time, A.time);
    end
    need_full(V.Vx, 'uv_slice', 'out.uv and out.modes');
    V.units = 'm s^{-1}';  V.short = 'u''';
    V.label = sprintf('baroclinic M_2 velocity at %g m', mopt.z);

  case 'uv_mean'
    % The observed total flow, averaged over the depths that reported in that
    % hour. See the header: the baroclinic depth mean is zero by construction.
    for i = 1:nm
        Mh = field_or(R{i}, 'Mh');
        if isempty(Mh) || ~isfield(Mh,'u'), continue; end
        V.Vx(i,:) = snap(depthmean(Mh.u, Mh.depth), R{i}.time, A.time);
        V.Vy(i,:) = snap(depthmean(Mh.v, Mh.depth), R{i}.time, A.time);
    end
    need_full(V.Vx, 'uv_mean', 'out.Mh');
    V.units = 'm s^{-1}';  V.short = 'u';
    V.label = 'observed velocity, depth mean over the sampled range';

  otherwise
    error('array_map_movie:vector', ...
        ['Unknown mopt.vector "%s". One of: flux, flux_mean, bt, uv_slice, ' ...
         'uv_mean.'], mopt.vector);
end

V.mag = hypot(V.Vx, V.Vy);
end

%% ======================================================================
function S = get_scalar(R, A, mopt)
% One scalar per mooring per frame, for the dot. Sampled at mopt.scalar_z
% (nearest reporting sensor) or averaged over sensors when that is empty.
known = {'T','S','chi','eps','N2','spd','dTdz'};
nm = numel(A.name);  nt = numel(A.time);
f  = char(mopt.scalar);
S  = struct('val', nan(nm,nt), 'name', f);

% chi and eps span decades, so they are log-scaled unless told otherwise. A
% linear colour scale on eps is one bright dot and four black ones.
if isempty(mopt.scalar_log)
    S.log = any(strcmpi(f, {'chi','eps','N2','chiSSF','epsSSF'}));
else
    S.log = logical(mopt.scalar_log);
end

got = false;
for i = 1:nm
    Mh = field_or(R{i}, 'Mh');
    if isempty(Mh) || ~isfield(Mh, f), continue; end
    got = true;
    X = Mh.(f);
    if isempty(mopt.scalar_z)
        x = mean(X, 1, 'omitnan');
    else
        % Nearest sensor by its record-mean depth, among those that report the
        % variable at all. Interpolating a scalar across 40 m of sensor spacing
        % would invent a profile this array does not have.
        zb = mean(Mh.depth, 2, 'omitnan');
        zb(all(~isfinite(X), 2)) = NaN;
        [~, j] = min(abs(zb - mopt.scalar_z));
        x = X(j,:);
        fprintf('  %s: %s taken from the sensor at %.0f m (asked %.0f m)\n', ...
                A.name{i}, f, zb(j), mopt.scalar_z);
    end
    if S.log, x(x <= 0) = NaN;  x = log10(x); end
    S.val(i,:) = snap(x, R{i}.time, A.time);
end

if ~got
    error('array_map_movie:noscalar', ...
        ['No mooring carries "%s". The scalar comes off out.Mh, which only ' ...
         'exists when the pipeline ran with opt.keep = ''full'' -- the array ' ...
         'loop in turner_fleat_IT_flux uses ''flux''. Known fields: %s.'], ...
         f, strjoin(known, ', '));
end

% Pooled percentiles, so one mooring cannot set the scale for the array. 2-98
% rather than min-max for the same reason the arrow scale uses the 99th.
v = S.val(isfinite(S.val));
if isempty(v)
    error('array_map_movie:scalarempty', '"%s" is finite nowhere.', f);
end
S.clim = [prctile(v,2) prctile(v,98)];
if diff(S.clim) <= 0, S.clim = S.clim(1) + [-0.5 0.5]; end

S.cmap  = pick_cmap(cmap_for(f), 256);
S.label = scalar_label(f, S.log, mopt.scalar_z);
end

%% ======================================================================
function B = get_bathy(R, mopt)
% Subset, decimate, and band. Returns the band INDEX as well as the depth,
% because the discrete colouring contours the index. See the header.
bf = mopt.bathy;
if isempty(bf)
    for i = 1:numel(R)
        c = field_or(R{i}, 'cfg');
        if ~isempty(c) && isfield(c,'bathy') && ~isempty(c.bathy)
            bf = c.bathy;  break
        end
    end
end
if isempty(bf) || ~isfile(bf)
    error('array_map_movie:nobathy', ...
        ['No bathymetry file. Pass mopt.bathy, or run with configs that set ' ...
         'cfg.bathy (fleat_config_base does). Tried: "%s"'], bf);
end

bl = double(ncread(bf,'lat'));
bo = double(ncread(bf,'lon'));
bz = double(ncread(bf,'z'));                 % [lon x lat], positive down

il = bl >= mopt.box(3) & bl <= mopt.box(4);
io = bo >= mopt.box(1) & bo <= mopt.box(2);
if ~any(il) || ~any(io)
    error('array_map_movie:box', ...
        'mopt.box is outside the bathymetry, which spans lon %.2f-%.2f, lat %.2f-%.2f.', ...
        min(bo), max(bo), min(bl), max(bl));
end
bl = bl(il);  bo = bo(io);  bz = bz(io, il);

d = max(1, round(numel(bo)/500));            % contour crawls on the full grid
bo = bo(1:d:end);  bl = bl(1:d:end);  bz = bz(1:d:end, 1:d:end);

B.lon = bo;  B.lat = bl;  B.z = bz.';
B.lev = mopt.levels;
if isempty(B.lev)
    % Tight over the reef and the shelf break, loose in the basin. The moorings
    % sit at 434, 880, 1515, 1666 and 3390 m, so every one of those falls in a
    % band of its own rather than being lumped with the basin floor.
    B.lev = [0 20 50 100 200 300 450 600 900 1200 1550 1700 2200 2800 3500 4500];
end
B.lev = unique(B.lev(:)).';
B.idx = discretize(B.z, B.lev);              % NaN outside, i.e. land and abyss

fprintf('bathymetry: %s\n  %d x %d after subsetting and decimating by %d, %.0f-%.0f m in box\n', ...
        bf, numel(bo), numel(bl), d, min(B.z(:)), max(B.z(:)));
fprintf('  %d discrete bands: %s m\n', numel(B.lev)-1, strjoin(compose('%g', B.lev), ' '));
if any(B.z(:) > B.lev(end))
    warning('array_map_movie:deeper', ...
        ['The box reaches %.0f m but the deepest band edge is %.0f m, so the ' ...
         'deepest water is uncoloured. Extend mopt.levels.'], ...
         max(B.z(:)), B.lev(end));
end
end

%% ======================================================================
function dot = place_dots(lon, lat, mdir, mopt)
% One pass, crowded moorings first, 24 candidate angles each. Score = distance
% to the nearest obstacle, in aspect-corrected degrees, with the arrow treated
% as an obstacle along its own direction. See LAYOUT in the header.
nm  = numel(lon);
c   = cosd(mean(mopt.box(3:4)));
X   = lon(:)*c;   Y = lat(:);                % the metric both are measured in
r   = mopt.dot_r;
ang = linspace(0, 2*pi, 25);  ang(end) = [];

% Crowded first: nearest-neighbour distance, ascending.
nn = inf(nm,1);
for i = 1:nm
    for j = 1:nm
        if j ~= i && isfinite(X(j)) && isfinite(X(i))
            nn(i) = min(nn(i), hypot(X(i)-X(j), Y(i)-Y(j)));
        end
    end
end
[~, order] = sort(nn, 'ascend');

dot = nan(nm,2);
P   = [X Y];                                  % obstacles: every mooring marker
bx  = mopt.box(1:2)*c;   by = mopt.box(3:4);  % the box, in the same metric
for i = order(:).'
    if ~isfinite(X(i)) || ~isfinite(Y(i)), continue; end
    % A circle of radius r in the ASPECT-CORRECTED metric, so it is a circle on
    % screen. r*c*cos here instead of r*cos would flatten it zonally by cos(lat)
    % and the candidate angles would no longer be equally spaced in distance.
    cx = X(i) + r*cos(ang);   cy = Y(i) + r*sin(ang);

    best = -inf;  bk = 1;
    for k = 1:numel(ang)
        % inside the box, with a margin of one dot radius
        if cx(k) < bx(1)+r || cx(k) > bx(2)-r || ...
           cy(k) < by(1)+r || cy(k) > by(2)-r
            continue
        end
        dmin = min(hypot(P(:,1)-cx(k), P(:,2)-cy(k)));
        % the arrow leaves its mooring along mdir(i); penalise sitting on it
        if isfinite(mdir(i))
            dth  = abs(angdiff_(ang(k), mdir(i)));
            dmin = min(dmin, r*dth/(pi/4));   % full credit beyond 45 deg off
        end
        if dmin > best, best = dmin;  bk = k; end
    end
    if ~isfinite(best)
        % Every candidate fell outside the box: the mooring is within one dot
        % radius of an edge. Leave it unplaced rather than drawing a dot off
        % the map, where manual axis limits would hide it anyway.
        warning('array_map_movie:nodot', ...
            ['No room for mooring %d''s scalar dot inside mopt.box at radius ' ...
             '%.3f deg. Shrink mopt.dot_r or widen the box.'], i, r);
        continue
    end
    dot(i,:) = [cx(bk)/c, cy(bk)];
    P(end+1,:) = [cx(bk) cy(bk)];             %#ok<AGROW> now an obstacle too
end

sep = inf;
for i = 1:nm
    for j = i+1:nm
        if all(isfinite(dot([i j],1)))
            sep = min(sep, hypot(diff(dot([i j],1))*c, diff(dot([i j],2))));
        end
    end
end
fprintf('scalar dots: radius %.3f deg, closest pair %.3f deg (%.1f km)\n', ...
        r, sep, sep*111.2);
end

%% ======================================================================
function rgb = val2rgb(v, cl, cmap)
% Values to RGB by hand, because the axes colormap belongs to the bathymetry.
% NaN comes back light grey: a dot with no data still marks its mooring.
v   = v(:);
n   = size(cmap,1);
f   = (v - cl(1)) / (cl(2) - cl(1));
f   = min(max(f, 0), 1);                       % clip, do not drop
k   = round(f*(n-1)) + 1;
rgb = repmat([.85 .85 .85], numel(v), 1);
ok  = isfinite(k);
rgb(ok,:) = cmap(k(ok),:);
end

%% ======================================================================
function y = snap(x, t, tq)
% Put one series on the common clock by rounding, never by interpolating --
% the same rule and the same reason as ALIGN_FLUX_RECORDS: these are already
% hourly means, and interpolation would invent values across gaps and ends.
y = nan(1, numel(tq));
if isempty(x) || isempty(t), return; end
x = x(:).';  t = t(:).';
if numel(tq) < 2, return; end
step = tq(2) - tq(1);
k    = round((t - tq(1))/step) + 1;
in   = k >= 1 & k <= numel(tq) & isfinite(k);
y(k(in)) = x(in);
end

%% ======================================================================
function m = depthmean(X, D)
% Mean over the sensors that reported, hour by hour. Trapezoid in depth would
% be better but needs a sorted, gap-free column; with one or two sensors
% dropping in and out the weights move more than the signal does.
m = mean(X, 1, 'omitnan');
m(all(~isfinite(X) | ~isfinite(D), 1)) = NaN;
end

%% ======================================================================
function need_full(Vx, what, fieldname)
if ~any(isfinite(Vx(:)))
    error('array_map_movie:needfull', ...
        ['mopt.vector = ''%s'' needs %s, which the pipeline only attaches ' ...
         'when it runs with opt.keep = ''full''. The array loop in ' ...
         'turner_fleat_IT_flux uses ''flux'', so rerun it with ''full'' (it ' ...
         'is much bigger in memory) or use mopt.vector = ''flux''.'], ...
         what, fieldname);
end
end

%% ======================================================================
function s = lbl(R, A, i)
s = A.name{i};
c = field_or(R{i}, 'cfg');
if ~isempty(c) && isfield(c,'unverified') && ~isempty(c.unverified)
    s = [s '*'];
end
end

%% ======================================================================
function s = unverified_note(R, A)
s = '';
for i = 1:numel(A.name)
    c = field_or(R{i}, 'cfg');
    if ~isempty(c) && isfield(c,'unverified') && ~isempty(c.unverified)
        s = [s sprintf('%s* unverified: %s.   ', A.name{i}, ...
                       strjoin(c.unverified, ', '))];  %#ok<AGROW>
    end
end
if ~isempty(s), s = []; end
end

%% ======================================================================
function c = cmap_for(f)
switch lower(f)
  case 't',          c = 'thermal';
  case 's',          c = 'haline';
  case {'chi','eps'},c = 'amp';
  case 'n2',         c = 'dense';
  otherwise,         c = 'amp';
end
end

%% ======================================================================
function s = scalar_label(f, islog, z)
switch lower(f)
  case 't',    s = 'temperature [\circC]';
  case 's',    s = 'salinity';
  case 'chi',  s = '\chi [K^2 s^{-1}]';
  case 'eps',  s = '\epsilon [W kg^{-1}]';
  case 'n2',   s = 'N^2 [s^{-2}]';
  otherwise,   s = f;
end
if islog, s = ['log_{10} ' s]; end
if isempty(z), s = [s '  (sensor mean)'];
else,          s = [s sprintf('  (nearest sensor to %g m)', z)];
end
end

%% ======================================================================
function cmap = pick_cmap(name, n)
% cmocean if it is on the path, a built-in if it is not. The figure should not
% fail to draw over a colormap.
try
    cmap = cmocean(name, n);
catch
    switch name
      case 'deep',    cmap = flipud(parula(n));
      case 'thermal', cmap = turbo(n);
      case 'haline',  cmap = parula(n);
      otherwise,      cmap = hot(n);
    end
end
end

%% ======================================================================
function v = round_nice(x)
% A reference arrow of 1, 2 or 5 times a power of ten, so the corner of the
% figure carries a number a reader can do arithmetic with.
%
% Guarded at the top because log10(0) is -Inf and log10(NaN) is NaN, either of
% which walks through the comparisons below and returns a zero-length reference
% arrow labelled "0". That happens when the caller supplies mopt.scale and the
% field is empty, so the arrow is the only thing left to notice it.
if ~isfinite(x) || x <= 0, v = 1; return; end
p = 10^floor(log10(x));
m = x/p;
if     m < 1.5, v = 1*p;
elseif m < 3.5, v = 2*p;
elseif m < 7.5, v = 5*p;
else,           v = 10*p;
end
end

%% ======================================================================
function d = angdiff_(a, b)
d = atan2(sin(a-b), cos(a-b));
end

%% ======================================================================
function v = field_or(s, f)
if isstruct(s) && isfield(s, f), v = s.(f); else, v = []; end
end

%% ======================================================================
function o = setdef(o, f, v)
if ~isfield(o, f), o.(f) = v; end
end

%% ======================================================================
function y = tern(c, a, b)
if c, y = a; else, y = b; end
end
