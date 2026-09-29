function h = psectime_mooring(structure, var1, var2, cont1, cont2, cont2h, ztype, zmax, cfg, tidx)
% PSECTIME_MOORING  Contour a mooring variable against depth or density vs time.
%
%   h = psectime_mooring(structure, var1, var2, cont1, cont2, cont2h, ztype, zmax, cfg)
%
%   structure  mooring struct (e.g. F2) with .time, .depth and the named fields
%   var1       fieldname contoured in colour
%   var2       fieldname contoured in black. [] to skip.
%   cont1      contour levels for var1. [] to autoset (~15 over robust range).
%   cont2      contour levels for var2. [] to autoset.
%   cont2h     heavy, labelled levels of var2. [] to autoset as every 5th of cont2.
%   ztype      'depth' | 'sigma'
%   zmax       max depth in m. [] for no cut. In 'sigma' mode this keeps density
%              surfaces whose MEAN depth is shallower than zmax, since depth is
%              not the axis there.
%   cfg        optional; cfg.mooring names the plot. Falls back to inputname.
%   tidx       optional; time indices to plot, e.g. 2654:3192. [] for full record.
%
%   Returns h with .ax, .cb, .c1, .c2 handles.
%
%   TWO TIME BASES. This struct carries 10-minute fields (time, depth, T, S, ...)
%   at 45072 columns and hourly fields (slope*) at 7512, matching slopetime. The
%   function picks the right time vector from each field's column count and, if
%   var1 and var2 differ, interpolates var2 onto var1's grid.
%
%   'SIGMA' MODE removes isopycnal heave. Each profile is Thorpe-sorted, then
%   var1 is interpolated onto a fixed set of target density surfaces, so a row of
%   the plot follows one isopycnal rather than one depth. The y-axis is sigma
%   increasing downward. Because contouring sigma against sigma is degenerate,
%   var2 is replaced by isopycnal DEPTH in this mode, which puts the heave back
%   as information rather than as distortion.

% author TEJ 9 September 2026, with help from Claude.

%% ---- arguments ------------------------------------------------------
if nargin < 5, cont2  = []; end
if nargin < 6, cont2h = []; end
if nargin < 7 || isempty(ztype), ztype = 'depth'; end
if nargin < 8, zmax = []; end
if nargin < 9 || isempty(cfg), cfg = struct(); end
if nargin < 10, tidx = []; end

if isfield(cfg,'mooring') && ~isempty(cfg.mooring)
    mname = cfg.mooring;
else
    mname = inputname(1);
    if isempty(mname), mname = 'mooring'; end
end
lat = 8.5;
if isfield(cfg,'lat'), lat = cfg.lat; end

ztype = lower(ztype);
if ~ismember(ztype, {'depth','sigma'})
    error('psectime_mooring:ztype','ztype must be ''depth'' or ''sigma'', got ''%s''.', ztype);
end

if ~isfield(structure, var1)
    error('psectime_mooring:novar1', ...
        '''%s'' is not a field of the mooring struct. Available: %s', ...
        var1, strjoin(fieldnames(structure)', ', '));
end

% sigma has to be supplied, not derived -- deriving it here would hide which
% pressure and which equation of state were used.
if strcmp(ztype,'sigma') && ~isfield(structure,'sigma')
    error('psectime_mooring:nosigma', ...
        ['ztype = ''sigma'' needs a .sigma field and this struct has none. Add it:\n' ...
         '    P = sw_pres(%s.depth, %.4g);\n' ...
         '    %s.sigma = sw_pden(%s.S, %s.T, P, 0) - 1000;'], ...
         mname, lat, mname, mname, mname);
end
if ~isempty(var2) && ~isfield(structure, var2) && ~strcmp(ztype,'sigma')
    error('psectime_mooring:novar2','''%s'' is not a field of the mooring struct.', var2);
end

%% ---- pick the matching time vector for each field -------------------
V1 = structure.(var1);
t1 = pick_time(structure, size(V1,2), var1);

if strcmp(ztype,'sigma')
    if ~isempty(var2) && ~strcmp(var2,'depth')
        fprintf('psectime_mooring: sigma mode -- replacing var2 (%s) with isopycnal depth.\n', var2);
    end
    var2 = 'depth';
end

V2 = []; t2 = [];
if ~isempty(var2)
    V2 = structure.(var2);
    t2 = pick_time(structure, size(V2,2), var2);
end

% put var2 on var1's time grid if they differ
if ~isempty(V2) && (numel(t1) ~= numel(t2) || any(abs(t1(:)-t2(:)) > 1e-9))
    V2 = interp1(t2(:), V2.', t1(:)).';
end

D = structure.depth;
tD = pick_time(structure, size(D,2), 'depth');
if numel(tD) ~= numel(t1) || any(abs(tD(:)-t1(:)) > 1e-9)
    D = interp1(tD(:), D.', t1(:)).';
end

%% ---- subset in time -------------------------------------------------
if ~isempty(tidx)
    tidx = tidx(tidx >= 1 & tidx <= numel(t1));
    if isempty(tidx)
        error('psectime_mooring:tidx','tidx is outside 1:%d.', numel(t1));
    end
    t1 = t1(tidx);  V1 = V1(:,tidx);  D = D(:,tidx);
    if ~isempty(V2), V2 = V2(:,tidx); end
end

%% ---- decimate in time, or contour will crawl ------------------------
% Contour tracing over ~45000 columns is minutes of work for a figure that
% cannot resolve them anyway. Aim for a few thousand columns.
ncol   = numel(t1);
target = 4000;
dec    = max(1, floor(ncol/target));
if dec > 1
    fprintf('psectime_mooring: decimating %d -> %d columns (every %d) for contouring.\n', ...
            ncol, numel(1:dec:ncol), dec);
    idx = 1:dec:ncol;
    t1  = t1(idx);  V1 = V1(:,idx);  D = D(:,idx);
    if ~isempty(V2), V2 = V2(:,idx); end
end

%% ---- build the vertical coordinate ----------------------------------
switch ztype
  case 'depth'
    [~, ord] = sort(mean(D, 2, 'omitnan'));
    D  = D(ord,:);  V1 = V1(ord,:);
    if ~isempty(V2), V2 = V2(ord,:); end

    Y = D;                                   % positive down
    X = repmat(t1(:).', size(Y,1), 1);
    ylab = '$z$ (m)';

    if ~isempty(zmax)
        keep = mean(D,2,'omitnan') <= zmax;
        Y = Y(keep,:); X = X(keep,:); V1 = V1(keep,:);
        if ~isempty(V2), V2 = V2(keep,:); end
    end

  case 'sigma'
    SIG = structure.sigma;
    tS  = pick_time(structure, size(SIG,2), 'sigma');
    if size(SIG,2) ~= numel(t1)
        SIG = interp1(tS(:), SIG.', t1(:)).';
    end

    % target surfaces from the record-mean profile, thinned to the sensor count
    smean = mean(SIG, 2, 'omitnan');
    smean = smean(isfinite(smean));
    sig_t = linspace(prctile(smean,2), prctile(smean,98), numel(smean)).';

    nt  = numel(t1);
    V1s = nan(numel(sig_t), nt);
    Zs  = nan(numel(sig_t), nt);
    for it = 1:nt
        s = SIG(:,it);  v = V1(:,it);  z = D(:,it);
        ok = isfinite(s) & isfinite(z);
        if nnz(ok) < 4, continue; end
        s = s(ok); v = v(ok); z = z(ok);
        [~,o] = sort(z);  s = s(o); v = v(o); z = z(o);
        ss = sort(s);                                   % Thorpe sort
        [su, iu] = unique(ss);
        if numel(su) < 2, continue; end
        V1s(:,it) = interp1(su, v(iu), sig_t);
        Zs(:,it)  = interp1(su, z(iu), sig_t);
    end

    V1 = V1s;  V2 = Zs;
    Y = repmat(sig_t, 1, nt);
    X = repmat(t1(:).', numel(sig_t), 1);
    ylab = '$\sigma_\theta$ [kg m$^{-3}$]';

    if ~isempty(zmax)
        keep = mean(Zs, 2, 'omitnan') <= zmax;
        Y = Y(keep,:); X = X(keep,:); V1 = V1(keep,:); V2 = V2(keep,:);
    end
end

%% ---- units, colormap, log handling ----------------------------------
[clab, cmapname, is_log, is_div] = var_style(var1);
if is_log
    V1 = log10(max(V1, realmin));
    clab = ['$\log_{10}$ ' clab];
end

%% ---- levels ---------------------------------------------------------
if isempty(cont1)
    lo = prctile(V1(:), 1);  hi = prctile(V1(:), 99);
    if is_div, m = max(abs([lo hi])); lo = -m; hi = m; end
    if ~(isfinite(lo) && isfinite(hi)) || lo == hi
        error('psectime_mooring:levels','Could not autoset cont1 for %s -- all NaN or constant.', var1);
    end
    cont1 = linspace(lo, hi, 15);
end
if ~isempty(V2) && isempty(cont2)
    lo = prctile(V2(:), 2);  hi = prctile(V2(:), 98);
    cont2 = linspace(lo, hi, 12);
end
if ~isempty(V2) && isempty(cont2h)
    cont2h = cont2(1:5:end);
end
bad = ~isfinite(V1);
V1  = min(max(V1, cont1(1)), cont1(end));
V1(bad) = NaN;
%% ---- draw -----------------------------------------------------------
h.ax = gca; hold(h.ax,'on'); box(h.ax,'on');

[~, h.c1] = contourf(h.ax, X, Y, V1, cont1, 'LineStyle','none');
colormap(h.ax, cmocean(cmapname, numel(cont1)-1));
caxis(h.ax, [min(cont1) max(cont1)]);

if ~isempty(V2)
    [~, h.c2] = contour(h.ax, X, Y, V2, cont2, 'k-', 'LineWidth', 0.4);
    if ~isempty(cont2h)
        [~, h.c2h] = contour(h.ax, X, Y, V2, cont2h, 'k-', 'LineWidth', 1.1, ...
                             'ShowText','on', 'LabelSpacing', 400);
    end
end

h.cb = colorbar(h.ax);
ylabel(h.cb, clab, 'Interpreter','latex', 'FontSize', 11);

ylabel(h.ax, ylab, 'Interpreter','latex', 'FontSize', 11);
if strcmp(ztype,'sigma')
    set(h.ax, 'YDir', 'reverse');       % denser at the bottom
end
set(h.ax, 'XLim', [t1(1) t1(end)]);
if t1(end) - t1(1) <= 60                          % short window: one tick per day
    set(h.ax, 'XTick', ceil(t1(1)):floor(t1(end)));
    datetick(h.ax, 'x', 'dd mmm', 'keeplimits', 'keepticks');
else                                              % full record: monthly
    datetick(h.ax, 'x', 'mmm yyyy', 'keeplimits');
end
set(h.ax, 'XTickLabelRotation', 45);

if isfield(cfg,'ylim') && ~isempty(cfg.ylim)
    set(h.ax, 'YLim', cfg.ylim);
elseif ~isempty(zmax)
    set(h.ax, 'YLim', [0 zmax]);
else
    set(h.ax, 'YLim', [0 max(Y(:))*1.02]);
end
set(h.ax, 'YDir', 'reverse');

title(h.ax, sprintf('%s:   %s - %s', mname, ...
      datestr(t1(1),'HH:MM dd mmm yyyy'), datestr(t1(end),'HH:MM dd mmm yyyy')), ...
      'FontWeight','normal', 'FontSize', 11);
end

% ======================================================================
function t = pick_time(S, n, name)
% Match a field's column count to the right time vector.
if isfield(S,'time') && numel(S.time) == n
    t = S.time(:).';
elseif isfield(S,'slopetime') && numel(S.slopetime) == n
    t = S.slopetime(:).';
else
    error('psectime_mooring:timebase', ...
        ['Field ''%s'' has %d columns, which matches neither time (%d) nor ' ...
         'slopetime (%d).'], name, n, ...
         numel(getfielddef(S,'time')), numel(getfielddef(S,'slopetime')));
end
end

function v = getfielddef(S, f)
if isfield(S,f), v = S.(f); else, v = []; end
end

% ======================================================================
function [lab, cmap, is_log, is_div] = var_style(v)
% Colorbar label (LaTeX, with units), cmocean palette, and whether the field
% wants a log transform or a symmetric diverging scale.
%
% Log for chi, eps and N2 because they span decades and a linear scale renders
% them as one flat colour. Diverging-and-symmetric for anything that changes
% sign, because an off-centre zero on a diverging map is actively misleading.
is_log = false;  is_div = false;

switch lower(v)
  case 't',            lab = '$T$ ($^\circ$C)';                        cmap = 'thermal';
  case 's',            lab = '$S$ (psu)';                              cmap = 'haline';
  case 'sigma',        lab = '$\sigma_\theta$ (kg m$^{-3}$)';          cmap = 'dense';
  case 'depth',        lab = '$z$ (m)';                                cmap = 'deep';
  case 'n2',           lab = '$N^2$ (s$^{-2}$)';                       cmap = 'tempo';  is_log = true;
  case 'dtdz',         lab = '$\partial T/\partial z$ ($^\circ$C m$^{-1}$)'; cmap = 'balance'; is_div = true;
  case 'spd',          lab = 'speed (m s$^{-1}$)';                     cmap = 'speed';
  case 'u',            lab = '$u$ (m s$^{-1}$)';                       cmap = 'balance'; is_div = true;
  case 'v',            lab = '$v$ (m s$^{-1}$)';                       cmap = 'balance'; is_div = true;
  case {'chi','chissf','chi_nomask'}
                       lab = '$\chi$ ($^\circ$C$^2$ s$^{-1}$)';        cmap = 'matter'; is_log = true;
  case {'eps','epsssf','eps_nomask'}
                       lab = '$\varepsilon$ (W kg$^{-1}$)';            cmap = 'matter'; is_log = true;
  case {'slope','slope3','slope20020','slope_corspec'}
                       lab = 'spectral slope';                         cmap = 'balance'; is_div = true;
  case 'eta',   lab = '$\eta$ (m, + up)';   cmap = 'balance'; is_div = true;
  case 'ziso',  lab = 'isopycnal depth (m)';  cmap = 'deep';
  otherwise
    lab  = strrep(v,'_','\_');
    cmap = 'thermal';
    warning('psectime_mooring:unknownvar', ...
        'No style entry for ''%s''; using a plain label and thermal.', v);
end
end