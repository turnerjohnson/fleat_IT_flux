function E = plot_tidal_ellipse(u, v, t, cfg, opt)
% PLOT_TIDAL_ELLIPSE  Tidal current ellipses by rotary decomposition.
%
%   E = plot_tidal_ellipse(u, v, t, cfg, opt)
%   E = plot_tidal_ellipse(F2_h, cfg, opt)        % struct shorthand
%
%   u, v  [n x nt]  velocity pair [m/s]. The ROWS can be anything: sensors at
%                   depth, modal amplitudes from fit_uv_modes, or a single
%                   depth-averaged series. NaNs allowed.
%   t     [1 x nt]  time in DAYS (datenum), as F2_h.time.
%   cfg             from f2_config -- uses .M2_cpd and .mooring. Not .lat: the
%                   fit is a harmonic at a given frequency and never touches
%                   the position, which is why plot_tidal_ellipses_latlon can
%                   hand it one config while looping over every mooring.
%
%   opt .freq_cpd   constituent to fit [cpd]. Default cfg.M2_cpd.
%       .rowkind    'mode' | 'depth' | 'row'. Sets labels and colours, and for
%                   'mode' it labels row 1 as mode 0 (barotropic), matching the
%                   uhat/vhat convention out of fit_uv_modes.
%       .z          [n x 1] depth of each row [m], for 'depth' labels.
%       .rows       subset of rows to use. Default all.
%       .segments   split the record into this many chunks and fit each one,
%                   to show how far the ellipse is from stationary. 0 = off.
%       .layout     'overlay' (default) | 'grid'
%       .minobs     min finite samples a row needs. Default 48 (2 days hourly).
%       .plot       false returns E without drawing anything, for callers that
%                   want the rotary fit but their own figure (see
%                   plot_tidal_ellipses_latlon). Default true.
%       .t0         phase reference epoch, in the same time base as t.
%                   Default 0, i.e. datenum 0, so phases are comparable
%                   between rows, moorings and deployments.
%
%   E.major, E.minor  [n x 1]  semi-axes [m/s]. MINOR IS SIGNED: positive means
%                              counterclockwise rotation, negative clockwise.
%   E.incl            [n x 1]  major-axis orientation [deg CCW from east), 0-180
%   E.pha             [n x 1]  phase [deg]: the value of omega*(t-t0) at which
%                              the current points along +incl and is fastest.
%   E.ellip           [n x 1]  minor/major, in [-1 1]. 0 rectilinear, +-1 circular.
%   E.U, E.V          [n x 1]  complex amplitudes, u = Re(U e^{i omega tau}).
%   E.varexp          [n x 1]  fraction of each row's variance the single
%                              harmonic explains. READ THIS -- see below.
%   E.seg             per-segment major axis, if opt.segments > 0.
%
%   ------------------------------------------------------------------
%   WHAT IS ACTUALLY BEING COMPUTED
%   ------------------------------------------------------------------
%   Write each component as a single harmonic at frequency omega,
%
%       u(t) = Re{U e^{i omega tau}},   v(t) = Re{V e^{i omega tau}}
%
%   with U, V complex and tau = t - t0. Then w = u + iv splits into two
%   counter-rotating circles (the rotary decomposition),
%
%       w(t) = W+ e^{i omega tau} + W- e^{-i omega tau}
%       W+ = (U + iV)/2,    W- = conj(U - iV)/2
%
%   and the tip of w traces an ellipse. The two circles line up once per half
%   cycle, giving semi-major = |W+| + |W-|, and oppose, giving semi-minor =
%   |W+| - |W-|. That difference is signed on purpose: |W+| > |W-| means the
%   counterclockwise circle wins and the current veers CCW. Orientation and
%   phase are the half-sum and half-difference of the two circles' arguments.
%
%   Verified against brute-force sampling of the hodograph: major reproduces
%   max|w|, |minor| reproduces min|w|, the sign of minor matches the sign of
%   the swept area, and at omega*tau = pha the vector lies on the +incl axis
%   at maximum speed.
%
%   ------------------------------------------------------------------
%   THIS IS NOT HARMONIC ANALYSIS, AND THE PHASE IS NOT A GREENWICH PHASE
%   ------------------------------------------------------------------
%   One frequency is fit by least squares. There is no nodal correction, no
%   equilibrium argument, and no attempt to separate M2 from its neighbours.
%   Over a six-month record M2 and S2 are well separated (Rayleigh needs 14.8
%   days) so the fit is clean, but the phase is referenced to opt.t0 in your
%   own time base, not to the astronomical argument. It is consistent across
%   rows and records, which is what you need to compare moorings; it is not a
%   harmonic constant you can publish as one. You have t_tide on your path --
%   use it if you need true constants with confidence intervals.
%
%   E.varexp is the check that matters. Feeding this RAW velocity gives a low
%   varexp because everything that is not M2 counts as residual; that is
%   expected and not a failure. Feeding it the already-bandpassed uhat_M2 from
%   fit_uv_modes should give varexp near 1, and if it does not, the amplitude
%   series is not behaving like a single narrowband harmonic -- which usually
%   means the fit it came from is unstable rather than that the tide is odd.
%
%   ------------------------------------------------------------------
%   STATIONARITY
%   ------------------------------------------------------------------
%   A single ellipse over the whole record assumes the tide is stationary. It
%   is not: spring-neap alone modulates the amplitude by a factor of ~2, and
%   the baroclinic tide is far less predictable than the barotropic one
%   because it depends on a stratification that moves. opt.segments fits each
%   chunk separately and draws the family faintly behind the mean ellipse. If
%   those scatter widely, the single ellipse is a summary statistic and not a
%   description, and the spread is the honest error bar to quote.

%% ---- arguments ------------------------------------------------------
if isstruct(u)
    S = u;
    if nargin >= 3, opt = t;   else, opt = struct(); end
    if nargin >= 2, cfg = v;   else, cfg = struct(); end
    if ~all(isfield(S, {'u','v','time'}))
        error('plot_tidal_ellipse:struct', ...
              'Struct input needs .u, .v and .time (as F2_h does).');
    end
    u = S.u;  v = S.v;  t = S.time;
    if ~isfield(opt,'rowkind'), opt.rowkind = 'depth'; end
    if ~isfield(opt,'z') && isfield(S,'depth')
        opt.z = mean(S.depth, 2, 'omitnan');
    end
else
    % nargin refers to the ARRAY form here; in the struct form cfg and opt
    % have already been shifted out of v and t above and must not be reset
    if nargin < 4, cfg = struct(); end
    if nargin < 5, opt = struct(); end
end
if isempty(cfg), cfg = struct(); end
if isempty(opt), opt = struct(); end

if ~isfield(cfg,'M2_cpd'),   cfg.M2_cpd = 1.9323;      end
if ~isfield(opt,'freq_cpd'), opt.freq_cpd = cfg.M2_cpd; end
if ~isfield(opt,'rowkind'),  opt.rowkind = 'row';       end
if ~isfield(opt,'rows'),     opt.rows    = [];          end
if ~isfield(opt,'segments'), opt.segments = 0;          end
if ~isfield(opt,'layout'),   opt.layout  = 'overlay';   end
if ~isfield(opt,'minobs'),   opt.minobs  = 48;          end
if ~isfield(opt,'t0'),       opt.t0      = 0;           end
if ~isfield(opt,'plot'),     opt.plot    = true;        end
if ~isfield(opt,'z'),        opt.z       = [];          end

t = t(:).';
if size(u,2) ~= numel(t) || ~isequal(size(u), size(v))
    error('plot_tidal_ellipse:size', ...
          'u and v must be [n x nt] and match t (%d columns).', numel(t));
end

if ~isempty(opt.rows)
    u = u(opt.rows,:);  v = v(opt.rows,:);
    if ~isempty(opt.z), opt.z = opt.z(opt.rows); end
end
n = size(u,1);

if ~isempty(opt.z) && numel(opt.z) ~= n
    error('plot_tidal_ellipse:zsize', ...
          'opt.z has %d entries but there are %d rows to label.', numel(opt.z), n);
end
if strcmpi(opt.rowkind,'depth') && isempty(opt.z)
    warning('plot_tidal_ellipse:noz', ...
            'rowkind is "depth" but no opt.z given -- labelling by row index.');
    opt.rowkind = 'row';
end

mname = 'mooring';
if isfield(cfg,'mooring') && ~isempty(cfg.mooring), mname = cfg.mooring; end

w = 2*pi*opt.freq_cpd;              % rad/day, since t is in days

%% ---- fit one harmonic per row ---------------------------------------
[E.U, E.V, E.varexp] = deal(nan(n,1));
E.ngood = zeros(n,1);

for i = 1:n
    [E.U(i), E.V(i), E.varexp(i), E.ngood(i)] = fit_one(u(i,:), v(i,:), t, w, opt);
end

[E.major, E.minor, E.incl, E.pha] = uv_to_ellipse(E.U, E.V);
E.ellip = E.minor ./ E.major;
E.freq_cpd = opt.freq_cpd;
E.t0 = opt.t0;
E.rowkind = opt.rowkind;
E.z = opt.z;

%% ---- optional: is it stationary? ------------------------------------
if opt.segments > 1
    edges = round(linspace(1, numel(t)+1, opt.segments+1));
    E.seg.major = nan(n, opt.segments);
    E.seg.minor = nan(n, opt.segments);
    E.seg.incl  = nan(n, opt.segments);
    E.seg.U     = nan(n, opt.segments);
    E.seg.V     = nan(n, opt.segments);
    E.seg.tmid  = nan(1, opt.segments);
    for k = 1:opt.segments
        kk = edges(k):edges(k+1)-1;
        E.seg.tmid(k) = mean(t(kk));
        for i = 1:n
            [Uk, Vk] = fit_one(u(i,kk), v(i,kk), t(kk), w, opt);
            E.seg.U(i,k) = Uk;  E.seg.V(i,k) = Vk;
        end
        [E.seg.major(:,k), E.seg.minor(:,k), E.seg.incl(:,k)] = ...
            uv_to_ellipse(E.seg.U(:,k), E.seg.V(:,k));
    end
    E.seg.major_std = std(E.seg.major, 0, 2, 'omitnan');
end

%% ---- report ---------------------------------------------------------
lbl = row_labels(n, opt);
fprintf('\n%s tidal ellipses at %.4f cpd (%.2f h), %d rows, %.0f days\n', ...
        mname, opt.freq_cpd, 24/opt.freq_cpd, n, t(end)-t(1));
fprintf('%-14s %8s %8s %7s %7s %7s %6s %5s\n', ...
        '', 'major', 'minor', 'ellip', 'incl', 'pha', 'varexp', 'rot');
for i = 1:n
    if ~isfinite(E.major(i))
        fprintf('%-14s %8s  (too few good samples: %d)\n', lbl{i}, '--', E.ngood(i));
        continue
    end
    if E.minor(i) >= 0, rot = 'CCW'; else, rot = 'CW'; end
    fprintf('%-14s %8.4f %8.4f %7.2f %7.1f %7.1f %6.2f %5s\n', ...
            lbl{i}, E.major(i), E.minor(i), E.ellip(i), ...
            E.incl(i), E.pha(i), E.varexp(i), rot);
end
if opt.segments > 1
    fprintf('spread of semi-major across %d segments (1 sd): %s m/s\n', ...
            opt.segments, mat2str(round(E.seg.major_std',4)));
end
fprintf(['incl is deg CCW from east; minor > 0 is counterclockwise; pha is ' ...
         'referenced to t0 = %g\n'], opt.t0);

%% ---- plot -----------------------------------------------------------
if ~opt.plot, return; end
col = row_colors(n, opt);
phi = linspace(0, 2*pi, 361);

switch lower(opt.layout)
  case 'overlay'
    figure; hold on
    h = gobjects(n,1);
    for i = 1:n
        if ~isfinite(E.U(i)), continue; end
        if opt.segments > 1
            for k = 1:opt.segments
                plot(real(E.seg.U(i,k)*exp(1i*phi)), real(E.seg.V(i,k)*exp(1i*phi)), ...
                     '-', 'Color', fade(col(i,:)), 'LineWidth',0.5);
            end
        end
        h(i) = draw_one(E.U(i), E.V(i), phi, col(i,:), 1.8);
    end
    axis equal; grid on; box on
    xlabel('u (m s^{-1})'); ylabel('v (m s^{-1})')
    title(sprintf('%s  %.4f cpd current ellipses', mname, opt.freq_cpd))
    ok = isgraphics(h);
    if nnz(ok) <= 12
        legend(h(ok), lbl(ok), 'Location','eastoutside','Box','off')
    else
        if strcmpi(opt.rowkind,'depth')
            colormap(gca, cmocean('deep')); cb = colorbar;
            clim([min(opt.z) max(opt.z)]); cb.Direction = 'reverse';
            ylabel(cb, 'depth (m)');
        end
    end
    axlim = max(abs(axis))*1.02;
    xlim([-axlim axlim]); ylim([-axlim axlim]);
    plot(xlim,[0 0],'k-','LineWidth',0.5); plot([0 0],ylim,'k-','LineWidth',0.5);

  case 'grid'
    nc = ceil(sqrt(n)); nr = ceil(n/nc);
    figure('Position',[100 100 220*nc 220*nr],'Color','w')
    lim = max(E.major(isfinite(E.major)))*1.15;
    if isempty(lim) || ~isfinite(lim) || lim == 0, lim = 1; end
    for i = 1:n
        subplot(nr, nc, i); hold on
        plot([-lim lim],[0 0],'-','Color',[.8 .8 .8]);
        plot([0 0],[-lim lim],'-','Color',[.8 .8 .8]);
        if isfinite(E.U(i))
            if opt.segments > 1
                for k = 1:opt.segments
                    plot(real(E.seg.U(i,k)*exp(1i*phi)), real(E.seg.V(i,k)*exp(1i*phi)), ...
                         '-', 'Color', fade(col(i,:)), 'LineWidth',0.5);
                end
            end
            draw_one(E.U(i), E.V(i), phi, col(i,:), 1.6);
        end
        axis equal; xlim([-lim lim]); ylim([-lim lim]); box on
        title(lbl{i}, 'FontWeight','normal')
        if i == 1, xlabel('u (m s^{-1})'); ylabel('v (m s^{-1})'); end
    end
    sgtitle(sprintf('%s  %.4f cpd current ellipses', mname, opt.freq_cpd))

  otherwise
    error('plot_tidal_ellipse:layout','Unknown layout "%s".', opt.layout);
end
end

%% ======================================================================
function [U, V, varexp, ngood] = fit_one(ui, vi, ti, w, opt)
% least squares against cos and sin at the one frequency
U = NaN; V = NaN; varexp = NaN;
ok = isfinite(ui) & isfinite(vi);
ngood = nnz(ok);
if ngood < opt.minobs, return; end

tau = ti(ok) - opt.t0;
A   = [cos(w*tau(:)), sin(w*tau(:))];
if rcond(A'*A) < 1e-12, return; end

uu = ui(ok).' - mean(ui(ok));         % the harmonic carries no mean
vv = vi(ok).' - mean(vi(ok));
cu = A\uu;   cv = A\vv;
U  = cu(1) - 1i*cu(2);                % u = Re(U e^{i w tau})
V  = cv(1) - 1i*cv(2);

r  = [uu; vv] - [A*cu; A*cv];
s  = [uu; vv];
varexp = 1 - sum(r.^2)/max(sum(s.^2), eps);
end

%% ======================================================================
function [major, minor, incl, pha] = uv_to_ellipse(U, V)
% rotary decomposition. minor is SIGNED: positive = counterclockwise.
Wp = (U + 1i*V)/2;
Wm = conj(U - 1i*V)/2;
ap = abs(Wp);  am = abs(Wm);

major = ap + am;
minor = ap - am;

th = (angle(Wp) + angle(Wm))/2;
g  = (angle(Wm) - angle(Wp))/2;

% put inclination in [0,pi) and shift the phase by the same branch so the
% pair stays consistent -- adding pi to both leaves Wp and Wm untouched
k    = floor(th/pi);
incl = (th - k*pi) * 180/pi;
pha  = mod((g - k*pi) * 180/pi, 360);
end

%% ======================================================================
function h = draw_one(U, V, phi, c, lw)
ue = real(U*exp(1i*phi));
ve = real(V*exp(1i*phi));
h  = plot(ue, ve, '-', 'Color', c, 'LineWidth', lw);

% major axis, and an arrowhead giving the sense of rotation
[~, im] = max(hypot(ue, ve));
plot([-ue(im) ue(im)], [-ve(im) ve(im)], '-', 'Color', fade(c), 'LineWidth', 0.8);

ia = mod(im + 45, numel(phi)) + 1;              % a point a bit further round
du = real(1i*U*exp(1i*phi(ia)));                % tangent = direction of travel
dv = real(1i*V*exp(1i*phi(ia)));
sc = 0.18*hypot(ue(im), ve(im))/max(hypot(du,dv), eps);
quiver(ue(ia), ve(ia), du*sc, dv*sc, 0, ...
       'Color', c, 'LineWidth', lw, 'MaxHeadSize', 6);
end

%% ======================================================================
function c = fade(c)
% lighten toward white. Done this way rather than with an alpha channel
% because RGBA on a Line 'Color' is only semi-documented.
c = c + (1 - c)*0.78;
end

%% ======================================================================
function lbl = row_labels(n, opt)
lbl = cell(n,1);
switch lower(opt.rowkind)
  case 'mode'
    for i = 1:n
        if i == 1, lbl{i} = 'mode 0 (BT)';
        else,      lbl{i} = sprintf('mode %d', i-1);
        end
    end
  case 'depth'
    for i = 1:n, lbl{i} = sprintf('%.0f m', opt.z(i)); end
  otherwise
    for i = 1:n, lbl{i} = sprintf('row %d', i); end
end
end

%% ======================================================================
function col = row_colors(n, opt)
if strcmpi(opt.rowkind,'depth') && n > 6
    col = cmocean('deep',n);
    [~, o] = sort(opt.z(:));      % colour by depth, not by row order
    col(o,:) = col;
elseif n <= 7
    col = lines(n);
else
    col = turbo(n);
end
end
