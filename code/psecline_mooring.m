function h = psecline_mooring(structure, var1, cont1, ztype, zmax, cfg, tidx)
% PSECLINE_MOORING  One coloured line per instrument, at its actual depth.
%
%   h = psecline_mooring(structure, var1, cont1,  ztype, zmax, cfg)
%
%   Same arguments as PSECTIME_MOORING. The difference is rendering: instead of
%   contouring a field, each sensor is drawn as its own line, positioned at its
%   INSTANTANEOUS depth and coloured by the value it measured.
%
%   A mooring is ~40 discrete time series, not a field. contourf and pcolor both
%   interpolate between adjacent sensors, so a 148 m gap (741 -> 889 m) renders
%   as smoothly as the 10 m spacing in the bottom cluster and looks equally well
%   resolved. Here the gaps are gaps. Dead sensors are visibly absent, sparse
%   regions look sparse, and knockdown shows up as the lines moving.
%
%   ztype = 'sigma' plots one line per TARGET ISOPYCNAL, with y = that surface's
%   depth and colour = var1 on it. So the wiggle is heave and the colour is the
%   property riding on the surface. That is a different and more useful picture
%   than forcing sigma onto the y-axis, which would just give horizontal lines.
%
%   NOTE: var_style and pick_time below are copies of the ones inside
%   PSECTIME_MOORING. Worth factoring out into their own files eventually.

%% ---- arguments ------------------------------------------------------
if nargin < 3, cont1  = []; end
if nargin < 4 || isempty(ztype), ztype = 'depth'; end
if nargin < 5, zmax = []; end
if nargin < 6 || isempty(cfg), cfg = struct(); end
if nargin < 7, tidx = []; end


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
    error('psecline_mooring:ztype','ztype must be ''depth'' or ''sigma''.');
end
if ~isfield(structure, var1)
    error('psecline_mooring:novar1', ...
        '''%s'' is not a field. Available: %s', var1, strjoin(fieldnames(structure)', ', '));
end
if strcmp(ztype,'sigma') && ~isfield(structure,'sigma')
    error('psecline_mooring:nosigma', ...
        ['ztype = ''sigma'' needs a .sigma field. Add it:\n' ...
         '    P = sw_pres(%s.depth, %.4g);\n' ...
         '    %s.sigma = sw_pden(%s.S, %s.T, P, 0) - 1000;'], mname, lat, mname, mname, mname);
end

%% ---- fields and time bases ------------------------------------------
V1 = structure.(var1);
t1 = pick_time(structure, size(V1,2), var1);



D  = structure.depth;
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
end

%% ---- decimate --------------------------------------------------------
ncol   = numel(t1);
target = 4000;
dec    = max(1, floor(ncol/target));
if dec > 1
    fprintf('psecline_mooring: decimating %d -> %d columns (every %d).\n', ...
            ncol, numel(1:dec:ncol), dec);
    k  = 1:dec:ncol;
    t1 = t1(k);  V1 = V1(:,k);  D = D(:,k);

end
nt = numel(t1);

%% ---- vertical coordinate ---------------------------------------------
switch ztype
  case 'depth'
    [~, ord] = sort(mean(D, 2, 'omitnan'));
    D = D(ord,:);  V1 = V1(ord,:);
    Y = D;
    ylab = '$z$ (m)';
    if ~isempty(zmax)
        keep = mean(D,2,'omitnan') <= zmax;
        Y = Y(keep,:); V1 = V1(keep,:);
    end

  case 'sigma'
    SIG = structure.sigma;
    tS  = pick_time(structure, size(SIG,2), 'sigma');
    if size(SIG,2) ~= ncol
        SIG = interp1(tS(:), SIG.', t1(:)).';
    else
        SIG = SIG(:, 1:dec:ncol);
    end

    smean = mean(SIG, 2, 'omitnan');
    smean = smean(isfinite(smean));
    sig_t = linspace(prctile(smean,2), prctile(smean,98), numel(smean)).';

    V1s = nan(numel(sig_t), nt);
    Zs  = nan(numel(sig_t), nt);
    for it = 1:nt
        s = SIG(:,it);  v = V1(:,it);  z = D(:,it);
        ok = isfinite(s) & isfinite(z);
        if nnz(ok) < 4, continue; end
        s = s(ok); v = v(ok); z = z(ok);
        [~,o] = sort(z);  s = s(o); v = v(o); z = z(o);
        ss = sort(s);                                % Thorpe sort
        [su, iu] = unique(ss);
        if numel(su) < 2, continue; end
        V1s(:,it) = interp1(su, v(iu), sig_t);
        Zs(:,it)  = interp1(su, z(iu), sig_t);
    end
    V1 = V1s;  Y = Zs; 
    ylab = '$z$ of isopycnal (m)';
    if ~isempty(zmax)
        keep = mean(Zs,2,'omitnan') <= zmax;
        Y = Y(keep,:); V1 = V1(keep,:);
    end
    fprintf('psecline_mooring: sigma mode -- one line per isopycnal, y = its depth.\n');
end

%% ---- style, log, levels ----------------------------------------------
[clab, cmapname, is_log, is_div] = var_style(var1);
if is_log
    V1 = log10(max(V1, realmin));
    clab = ['$\log_{10}$ ' clab];
end

if isempty(cont1)
    lo = prctile(V1(:), 1);  hi = prctile(V1(:), 99);
    if is_div, m = max(abs([lo hi])); lo = -m; hi = m; end
    if ~(isfinite(lo) && isfinite(hi)) || lo == hi
        error('psecline_mooring:levels','Could not autoset cont1 for %s.', var1);
    end
    cont1 = linspace(lo, hi, 15);
end
bad = ~isfinite(V1);
V1  = min(max(V1, cont1(1)), cont1(end));
V1(bad) = NaN;

%% ---- draw -------------------------------------------------------------
h.ax = gca; hold(h.ax,'on'); box(h.ax,'on');

% Coloured lines. MATLAB has no native per-vertex line colour, so each row is
% drawn as a degenerate surface with EdgeColor interp -- the standard trick.
X = t1(:).';
for i = 1:size(V1,1)
    y = Y(i,:);  c = V1(i,:);
    if all(~isfinite(c)), continue; end          % dead sensor: leave it absent
    if range(c(isfinite(c))) == 0
        warning('psecline_mooring:flatrow', ...
            'Row %d is constant — check for a clamp applied over NaNs.', i);
    end
    surface(h.ax, [X;X], [y;y], zeros(2,nt), [c;c], ...
            'FaceColor','none', 'EdgeColor','interp', 'LineWidth', 1.1);
end

colormap(h.ax, cmocean(cmapname, numel(cont1)-1));
caxis(h.ax, [min(cont1) max(cont1)]);



h.cb = colorbar(h.ax);
ylabel(h.cb, clab, 'Interpreter','latex', 'FontSize', 11);
ylabel(h.ax, ylab, 'Interpreter','latex', 'FontSize', 11);
set(h.ax, 'XLim', [t1(1) t1(end)]);
if t1(end) - t1(1) <= 60                          % short window: one tick per day
    set(h.ax, 'XTick', ceil(t1(1)):floor(t1(end)));
    datetick(h.ax, 'x', 'dd mmm', 'keeplimits', 'keepticks');
else                                              % full record: monthly
    datetick(h.ax, 'x', 'mmm yyyy', 'keeplimits');
end
set(h.ax, 'XTickLabelRotation', 45);

nlive = nnz(any(isfinite(V1), 2));
if isfield(cfg,'ylim') && ~isempty(cfg.ylim)
    set(h.ax, 'YLim', cfg.ylim);
elseif ~isempty(zmax)
    set(h.ax, 'YLim', [0 zmax]);
else
    set(h.ax, 'YLim', [0 max(Y(:))*1.02]);
end
set(h.ax, 'YDir', 'reverse');
title(h.ax, sprintf('%s:   %s - %s   (%d of %d lines have data)', mname, ...
      datestr(X(1),'HH:MM dd mmm yyyy'), datestr(X(end),'HH:MM dd mmm yyyy'), ...
      nlive, size(V1,1)), 'FontWeight','normal', 'FontSize', 11);

end

% ======================================================================
function t = pick_time(S, n, name)
if isfield(S,'time') && numel(S.time) == n
    t = S.time(:).';
elseif isfield(S,'slopetime') && numel(S.slopetime) == n
    t = S.slopetime(:).';
else
    error('psecline_mooring:timebase', ...
        'Field ''%s'' has %d columns, matching neither time nor slopetime.', name, n);
end
end

% ======================================================================
function [lab, cmap, is_log, is_div] = var_style(v)
is_log = false;  is_div = false;
switch lower(v)
  case 't',     lab = '$T$ ($^\circ$C)';                cmap = 'thermal';
  case 's',     lab = '$S$ (psu)';                      cmap = 'haline';
  case 'sigma', lab = '$\sigma_\theta$ (kg m$^{-3}$)';  cmap = 'dense';
  case 'depth', lab = '$z$ (m)';                        cmap = 'deep';
  case 'n2',    lab = '$N^2$ (s$^{-2}$)';               cmap = 'tempo';   is_log = true;
  case 'dtdz',  lab = '$\partial T/\partial z$ ($^\circ$C m$^{-1}$)'; cmap='balance'; is_div=true;
  case 'spd',   lab = 'speed (m s$^{-1}$)';             cmap = 'speed';
  case 'u',     lab = '$u$ (m s$^{-1}$)';               cmap = 'balance'; is_div = true;
  case 'v',     lab = '$v$ (m s$^{-1}$)';               cmap = 'balance'; is_div = true;
  case {'chi','chissf','chi_nomask'}
                lab = '$\chi$ ($^\circ$C$^2$ s$^{-1}$)'; cmap = 'matter'; is_log = true;
  case {'eps','epsssf','eps_nomask'}
                lab = '$\varepsilon$ (W kg$^{-1}$)';    cmap = 'matter';  is_log = true;
  case {'slope','slope3','slope20020','slope_corspec'}
                lab = 'spectral slope';                 cmap = 'balance'; is_div = true;
    case 'eta',   lab = '$\eta$ (m)';             cmap = 'balance'; is_div = true;
    otherwise
    lab = strrep(v,'_','\_');  cmap = 'thermal';
    warning('psecline_mooring:unknownvar','No style entry for ''%s''.', v);
end
end