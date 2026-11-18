function flux_movie(F, th, cfg, bathyfile, opt)
% FLUX_MOVIE  Flux vector on a bathymetry map, and a stick plot of the record.
%
%   flux_movie(F, th, cfg)                  % position and bathy from cfg
%   flux_movie(F, th, cfg, bathyfile, opt)
%
%   F          struct from compute_tidal_flux
%   th         [1 x nt] time, datenum
%   cfg        from f2_config or f4_config -- uses .lat, .lon, .H, .mooring
%              and .bathy. This is the whole reason the position lives in the
%              config: the call used to carry 8.5027 and 134.5918 as literals,
%              and the F4 driver inherited them from the F2 driver, so the F4
%              movie drew its arrow at F2 and contoured F2's 1515 m isobath
%              around an 880 m mooring.
%   bathyfile  netCDF with lat, lon, z (z is [lon x lat]). Default cfg.bathy.
%   opt .frames    time indices to animate (default: one per day)
%       .outfile   default 'flux_<cfg.mooring>.mp4'. '' animates without writing.
%       .fps       default 10
%       .box       [lonmin lonmax latmin latmax] map extent
%       .scale     degrees of arrow per kW/m (default: auto from the data)
%       .isobath   contour drawn in white through the mooring. Default cfg.H.
%       .stick     true (default) also draws the stick plot in a second figure
%
%   ONE MOORING, ONE ARROW. A lat/lon field of flux needs many moorings --
%   Alford's map has 60. With F2 alone you get one vector that rotates and
%   changes length through the record. Process F1, F3, F4, F6 the same way and
%   the map fills in.
%
%   The STICK PLOT is usually the more readable of the two: the whole record in
%   one static figure, with each stick showing direction and magnitude at that
%   time. An animation of a single arrow is harder to read and harder to put in
%   a paper. Both are produced here.

%% ---- arguments --------------------------------------------------------
% The old form passed lat0 and lon0 as arguments 3 and 4. That is a hard error
% rather than a shim: remapping it silently would read argument 4 as the
% bathymetry file, and a six-argument call cannot reach this line at all.
if ~isstruct(cfg)
    error('flux_movie:oldform', ...
        ['flux_movie now takes the mooring config where lat0 was:\n' ...
         '    flux_movie(F, th, cfg, bathyfile, opt)\n' ...
         'cfg.lat, cfg.lon, cfg.H and cfg.bathy come from f2_config or ' ...
         'f4_config, so the position is no longer written at the call site.']);
end

if nargin < 5 || isempty(opt), opt = struct(); end
if nargin < 4 || isempty(bathyfile)
    if ~isfield(cfg,'bathy') || isempty(cfg.bathy)
        error('flux_movie:nobathy', ...
            'No bathymetry: pass one, or set cfg.bathy in the mooring config.');
    end
    bathyfile = cfg.bathy;
end

for f = {'lat','lon'}
    if ~isfield(cfg, f{1}) || ~isscalar(cfg.(f{1})) || ~isfinite(cfg.(f{1}))
        error('flux_movie:nopos', ...
            'cfg.%s must be a finite scalar -- see the site block in f2_config.', f{1});
    end
end
lat0 = cfg.lat;  lon0 = cfg.lon;
if ~isfield(cfg,'mooring'), cfg.mooring = ''; end

if ~isfield(opt,'outfile')
    if isempty(cfg.mooring)
        opt.outfile = 'flux.mp4';
    else
        opt.outfile = sprintf('flux_%s.mp4', cfg.mooring);
    end
end
if ~isfield(opt,'fps'),     opt.fps     = 10;   end
if ~isfield(opt,'stick'),   opt.stick   = true; end
if ~isfield(opt,'box'),     opt.box = [lon0-1.2 lon0+1.2 lat0-1.2 lat0+1.2]; end
if ~isfield(opt,'isobath')
    if isfield(cfg,'H'), opt.isobath = cfg.H; else, opt.isobath = []; end
end

if isempty(cfg.mooring), name_prefix = ''; else, name_prefix = [cfg.mooring '  ']; end

%% ---- bathymetry -------------------------------------------------------
bl = double(ncread(bathyfile,'lat'));
bo = double(ncread(bathyfile,'lon'));
bz = double(ncread(bathyfile,'z'));          % [lon x lat]

il = bl >= opt.box(3) & bl <= opt.box(4);
io = bo >= opt.box(1) & bo <= opt.box(2);
if ~any(il) || ~any(io)
    error('flux_movie:box', ...
        'Box is outside the bathymetry, which spans lon %.2f-%.2f, lat %.2f-%.2f.', ...
        min(bo), max(bo), min(bl), max(bl));
end
bl = bl(il);  bo = bo(io);  bz = bz(io, il);

% decimate -- 1991x1991 is 4M points and contour will crawl
d = max(1, round(numel(bo)/400));
bo = bo(1:d:end);  bl = bl(1:d:end);  bz = bz(1:d:end, 1:d:end);
fprintf('bathymetry: %d x %d after subsetting and decimating by %d\n', ...
        numel(bo), numel(bl), d);

BZ = bz.';


%% ---- arrow scaling ----------------------------------------------------
mag_kw = F.mag/1e3;
if ~isfield(opt,'scale') || isempty(opt.scale)
    ref = prctile(mag_kw, 95);
    if ~isfinite(ref) || ref == 0, ref = 1; end
    opt.scale = 0.5/ref;                     % 95th pct arrow spans 0.5 deg
end
fprintf('arrow scale: %.3f deg per kW/m (95th pct flux = %.2f kW/m)\n', ...
        opt.scale, prctile(mag_kw,95));

if ~isfield(opt,'frames') || isempty(opt.frames)
    ok = find(isfinite(F.Fx) & isfinite(F.Fy));
    step = round(24/ (th(2)-th(1)) / 24);     % ~daily
    step = max(step, round(24));
    opt.frames = ok(1:step:end);
end
fprintf('%d frames\n', numel(opt.frames));

%% ---- movie ------------------------------------------------------------
fh = figure('Position',[80 80 760 700],'Color','w');
% Locked geometry: a figure resized mid-render changes the captured frame
% size, which shears the movie. See capture_frame.
set(fh, 'Resize', 'off', 'PaperPositionMode', 'auto');
FRAME_DPI = 150;   % capture resolution, dots per inch
fsz = [];          % frame size, locked by the first capture_frame call
ax = axes(fh); hold(ax,'on'); box(ax,'on');

contourf(ax, bo, bl, BZ, 0:250:6000, 'LineStyle','none'); hold on;
colormap(ax, cmocean('deep'));       % flip: shallow light, deep dark
cb = colorbar(ax);
ylabel(cb, 'depth (m)');
caxis(ax, [0 6000]);
set(cb, 'Direction', 'reverse');             % 0 at the top, 6000 at the bottom
daspect(ax, [1 cosd(lat0) 1]);
xlabel(ax,'longitude (\circE)'); ylabel(ax,'latitude (\circN)');
% The mooring's own isobath, so the arrow can be read against the depth it
% actually sits on: 1515 m at F2, 880 m at F4.
if ~isempty(opt.isobath) && isfinite(opt.isobath)
    contour(ax, bo, bl, BZ, opt.isobath*[1 1], 'w-', 'LineWidth', 1);
end

% reference arrow, fixed
xr = opt.box(1) + 0.08*diff(opt.box(1:2));
yr = opt.box(3) + 0.08*diff(opt.box(3:4));
quiver(ax, xr, yr, opt.scale/cosd(lat0), 0, 0, 'k', 'LineWidth',1.5, 'MaxHeadSize',2);
text(ax, xr, yr-0.06*diff(opt.box(3:4)), '1 kW m^{-1}', 'FontSize',9);

hq = quiver(ax, lon0, lat0, 0, 0, 0, 'Color',[0.85 0.1 0.1], ...
            'LineWidth',2.5, 'MaxHeadSize',1.5);
ht = title(ax, '', 'FontWeight','normal');

if ~isempty(opt.outfile)
    vw = VideoWriter(opt.outfile,'MPEG-4'); vw.FrameRate = opt.fps; open(vw);
end

for k = 1:numel(opt.frames)
    it = opt.frames(k);
    % degrees per kW/m, with the zonal component stretched for the projection
    du = F.Fx(it)/1e3 * opt.scale / cosd(lat0);
    dv = F.Fy(it)/1e3 * opt.scale;
    set(hq, 'UData', du, 'VData', dv);
    set(ht, 'String', sprintf('%s%s      |F| = %.2f kW m^{-1}, %.0f\\circ', ...
        name_prefix, datestr(th(it),'dd mmm yyyy'), mag_kw(it), F.dir(it)));
    drawnow;
    if ~isempty(opt.outfile)
        [im, fsz] = capture_frame(fh, FRAME_DPI, fsz);
        writeVideo(vw, im);
    end
end

if ~isempty(opt.outfile)
    close(vw);
    fprintf('wrote %s\n', opt.outfile);
end

%% ---- stick plot -------------------------------------------------------
% Convention: each stick starts on the zero line at its own time. The vertical
% component is northward flux in kW/m; the horizontal tilt encodes eastward
% flux, scaled so that the largest flux tilts by ~1/30 of the record length.
% Horizontal distance on this axis is time, so the tilt is a visual cue for
% direction rather than a quantity you read off the x-axis.
if opt.stick
    s  = opt.frames;
    xs = th(s);
    tilt_scale = (range(xs)/30) / max(mag_kw(s), [], 'omitnan');

    figure('Position',[80 80 1100 420],'Color','w'); hold on; box on
    quiver(xs, zeros(size(xs)), (F.Fx(s)/1e3)*tilt_scale, F.Fy(s)/1e3, 0, ...
           'Color',[0.2 0.3 0.6], 'LineWidth',0.8, 'ShowArrowHead','off');
    yline(0, 'k-');
    datetick('x','mmm','keeplimits');
    set(gca, 'XLim', [min(xs) max(xs)]);
    ylabel('flux (kW/m)');
    title([name_prefix 'stick plot: vertical = northward, tilt = eastward'], ...
          'FontWeight','normal');
end
end