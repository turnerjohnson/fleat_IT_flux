function B = tpxo_m2_bt(lon, lat, tvec, H, opt)
%TPXO_M2_BT  Barotropic M2 current at a point, from the TPXO10 atlas.
%
%   B = tpxo_m2_bt(lon, lat, tvec, H)
%   B = tpxo_m2_bt(lon, lat, tvec, H, opt)
%
%   lon, lat  the site, degrees. lon may be negative or 0-360.
%   tvec      datenums to predict at (hourly is fine, any spacing works)
%   H         water depth [m] used to turn TRANSPORT into VELOCITY. See
%             WHY THIS NEEDS A DEPTH -- it is the biggest error in the answer.
%   opt .file   the atlas .nc. Default ../../../data/u_m2_tpxo10_atlas_30_v2.nc
%       .npad   half-width of the node patch searched for water. Default 4.
%
%   B .u .v         [1 x nt]  predicted M2 barotropic velocity [m/s], east/north
%     .amp_u .amp_v scalar    amplitude [m/s]
%     .pha_u .pha_v scalar    Greenwich phase lag [deg]
%     .U .V         scalar    transport amplitude [m^2/s], what the file holds
%     .lon_u .lat_u .lon_v .lat_v   the nodes actually used
%     .moved_km     how far the nearest WATER node was from the site
%     .H .file .freq
%
%   ---------------------------------------------------------------------
%   WHAT THE FILE ACTUALLY CONTAINS
%   ---------------------------------------------------------------------
%   u_m2_tpxo10_atlas_30_v2.nc is an OTIS TRANSPORT file: uRe, uIm, vRe, vIm
%   are the complex amplitude of depth-integrated transport in cm^2/s, stored
%   as int32 with NO scale factor (checked against the nodes at F2 and F3).
%   So the velocity is
%
%       amp = |uRe + i*uIm| * 1e-4 / H          [m/s]
%
%   where 1e-4 converts cm^2/s to m^2/s and H is the water depth. The phase
%   convention is the file's own: GMT phase = atan2(-uIm, uRe), which is a
%   Greenwich phase LAG in the same sense t_tide uses, so it can be handed
%   straight to t_predic.
%
%   THE GRID IS STAGGERED. u lives on (lon_u, lat_u) and v on (lon_v, lat_v),
%   offset half a cell -- 1/60 deg, about 1.8 km here. Each component is read
%   on its own grid. Sampling v at the u node is a real error at a site on a
%   slope as steep as Palau's, and it is silent.
%
%   LAND IS ZERO, not NaN. The nearest node to a mooring on the island flank
%   can easily be dry, so a patch of +-opt.npad nodes is read and the nearest
%   node with nonzero transport wins. B.moved_km says how far that was; if it
%   is more than a few km, the number is not about your site any more.
%
%   ---------------------------------------------------------------------
%   WHY THIS NEEDS A DEPTH, AND WHY THAT IS THE WEAK POINT
%   ---------------------------------------------------------------------
%   The atlas stores transport, so velocity needs a depth, and the depth that
%   belongs in the denominator is TPXO'S OWN h at that node -- the h its
%   solution used. That lives in grid_tpxo10_atlas_30_v2.nc, which is not in
%   data/. Until it is, H is the mooring's recorded depth (cfg.H), and the
%   velocity is wrong by exactly the ratio of the two depths. At 1/30 deg on
%   the Velasco flank that ratio is not close to 1, and it is a SCALE error:
%   it stretches the amplitude without touching the phase or the shape, so
%   the time series looks perfectly healthy while being 20% or 50% off.
%
%   Download the grid file and pass its h to fix it properly.
%
%   ---------------------------------------------------------------------
%   AND IT IS ONE CONSTITUENT
%   ---------------------------------------------------------------------
%   This file is M2 alone (the 'con' variable reads 'm2'). A single
%   constituent CANNOT show a spring/neap cycle: the fortnight is the beat
%   between M2 and S2, so a pure M2 prediction has a constant envelope apart
%   from the 18.6-year nodal factor. If the question is whether the baroclinic
%   flux is modulated at the beat, this file can supply the M2 reference but
%   NOT the barotropic envelope to compare it against. For that, fetch at
%   least u_s2_tpxo10_atlas_30_v2.nc and sum the two predictions.
%
%   Nodal corrections (f and V+u) come from t_tide's t_predic via 'latitude',
%   evaluated at the centre of tvec, which is what OTPS does internally.
%
% author TEJ, with help from Claude
% 12 November 2026

if nargin < 5, opt = struct(); end
if ~isfield(opt,'npad') || isempty(opt.npad), opt.npad = 4; end
if ~isfield(opt,'file') || isempty(opt.file)
    % .../fleat_IT_flux/src/core -> .../research_fleat, then data/
    up = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    opt.file = fullfile(up, 'data', 'u_m2_tpxo10_atlas_30_v2.nc');
    if ~isfile(opt.file)
        alt = fullfile(up, 'fleat_IT_flux', 'data', 'u_m2_tpxo10_atlas_30_v2.nc');
        if isfile(alt), opt.file = alt; end
    end
end
if ~isfile(opt.file)
    error('tpxo_m2_bt:nofile', ...
        ['No atlas at %s. Pass opt.file, or put ' ...
         'u_m2_tpxo10_atlas_30_v2.nc in research_fleat/data.'], opt.file);
end
if isempty(H) || ~isfinite(H) || H <= 0
    error('tpxo_m2_bt:nodepth', ...
        ['H must be a positive depth in metres: the atlas holds TRANSPORT ' ...
         'and there is no velocity without it.']);
end

tvec = tvec(:)';
B.file = opt.file;  B.H = H;

%% ---- the two grids ---------------------------------------------------
lonq = mod(lon, 360);              % the atlas runs 0-360
lonu = ncread(opt.file, 'lon_u');  latu = ncread(opt.file, 'lat_u');
lonv = ncread(opt.file, 'lon_v');  latv = ncread(opt.file, 'lat_v');

[B.U, B.pha_u, B.lon_u, B.lat_u, du] = ...
    node(opt.file, 'uRe', 'uIm', lonu, latu, lonq, lat, opt.npad);
[B.V, B.pha_v, B.lon_v, B.lat_v, dv] = ...
    node(opt.file, 'vRe', 'vIm', lonv, latv, lonq, lat, opt.npad);
B.moved_km = max(du, dv);

%% ---- transport [cm^2/s] -> velocity [m/s] ---------------------------
B.amp_u = B.U * 1e-4 / H;
B.amp_v = B.V * 1e-4 / H;

%% ---- predict ---------------------------------------------------------
% freq in cycles per hour, from t_tide's own table so this file and t_predic
% cannot disagree about what M2 is.
B.freq = 0.0805114;
try
    c  = t_getconsts;
    ix = find(strcmp(cellstr(c.name), 'M2'), 1);
    if ~isempty(ix), B.freq = c.freq(ix); end
catch
    warning('tpxo_m2_bt:noconsts', ...
        't_getconsts unavailable; using M2 = %.7f cph.', B.freq);
end

B.u = t_predic(tvec, 'M2  ', B.freq, [B.amp_u 0 B.pha_u 0], 'latitude', lat);
B.v = t_predic(tvec, 'M2  ', B.freq, [B.amp_v 0 B.pha_v 0], 'latitude', lat);
B.u = real(B.u(:)');  B.v = real(B.v(:)');

fprintf(['tpxo M2 at %.4f, %.4f (H = %.0f m): |u| %.2f cm/s (G %.0f deg), ' ...
         '|v| %.2f cm/s (G %.0f deg)\n'], ...
        lon, lat, H, 100*B.amp_u, B.pha_u, 100*B.amp_v, B.pha_v);
if B.moved_km > 3
    fprintf(2, '  nearest WATER node is %.1f km away -- check this is still your site\n', ...
            B.moved_km);
end
end

% =======================================================================
function [amp, pha, lonn, latn, dkm] = node(f, vre, vim, lonv, latv, lon, lat, npad)
%NODE  Nearest water node of one component, and its amplitude and GMT phase.
%
%   Reads a (2*npad+1)^2 patch rather than the array: the file is 933 MB and
%   each of its four fields is 10800 x 5401.

[~, i0] = min(abs(lonv - lon));
[~, j0] = min(abs(latv - lat));
nx = numel(lonv);  ny = numel(latv);

ia = max(1, i0-npad);  ib = min(nx, i0+npad);
ja = max(1, j0-npad);  jb = min(ny, j0+npad);

Re = double(ncread(f, vre, [ia ja], [ib-ia+1, jb-ja+1]));
Im = double(ncread(f, vim, [ia ja], [ib-ia+1, jb-ja+1]));

% Land is stored as 0 in BOTH parts, so a node is water when either is nonzero.
wet = (Re ~= 0) | (Im ~= 0);
if ~any(wet(:))
    error('tpxo_m2_bt:dry', ...
        ['Every node within %d of (%.3f, %.3f) is land. Raise opt.npad if the ' ...
         'site really is this close to the coast.'], npad, lon, lat);
end

[LO, LA] = ndgrid(lonv(ia:ib), latv(ja:jb));
% Local flat-earth distance is plenty at 1/30 deg and avoids a dependency.
dx = (LO - lon) .* cosd(lat) * 111.195;
dy = (LA - lat) * 111.195;
d  = hypot(dx, dy);
d(~wet) = Inf;

[dkm, k] = min(d(:));
Z    = complex(Re(k), Im(k));
amp  = abs(Z);                                  % cm^2/s, transport
pha  = mod(atan2(-imag(Z), real(Z))*180/pi, 360);   % the file's GMT phase
lonn = LO(k);  latn = LA(k);
end
