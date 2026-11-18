function p = psectime_glance_mooring(F_h, cfg, temp_cont, sal_cont, uv_cont, time_lim, zlim, figbase)
% PSECTIME_GLANCE_MOORING  The four-panel look at a mooring, plus knockdown.
%
%   p = psectime_glance_mooring(F_h, cfg)
%   p = psectime_glance_mooring(F_h, cfg, temp_cont, sal_cont, uv_cont, time_lim)
%   p = psectime_glance_mooring(..., time_lim, zlim, figbase)
%
%   F_h        hourly mooring struct, after compute_density_mooring, so it has
%              .time .depth .T .S .u .v .sigma
%   cfg        mooring config. Only cfg.mooring and cfg.lat are read, and only
%              to title the panels; [] or a bare struct works.
%
%   temp_cont  contour levels for T.    [] autosets (~15 over a robust range)
%   sal_cont   contour levels for S.    []
%   uv_cont    contour levels for u AND v, so the two share a scale.  []
%   time_lim   time INDEX range, e.g. 2900:2999.  [] for the full record.
%   zlim       max depth [m], a SCALAR.  [] for the full column. Default [].
%   figbase    added to both figure numbers. Default 0, giving figures 1 and 2.
%
%   FIGURE figbase+1  T, S, u, v contoured against DEPTH, four stacked panels,
%                     with sigma contoured over the top in black.
%   FIGURE figbase+2  isopycnal depths, one line per target surface, over the
%                     WHOLE record whatever time_lim says -- that panel exists
%                     to show knockdown, and knockdown is episodic (Jul-Aug and
%                     Dec-Jan at F2), so a 100-hour window is the wrong view of
%                     it. zlim does apply, since cutting the depth axis does not
%                     hide when the array moved.
%
%   Returns p with the handle struct from each call: p.fig1, p.fig2, p.T, p.S,
%   p.u, p.v, p.sigma.
%
%   ---------------------------------------------------------------------
%   DEPTH ON THE AXIS, DENSITY ON TOP IN BLACK
%   ---------------------------------------------------------------------
%   These calls pass ztype = 'depth' and var2 = 'sigma', so the y-axis is depth
%   and the isopycnals are drawn over the colour as contours. Heave is therefore
%   visible as the black contours moving, which at a mooring knocked down 56 m
%   in the median and 63 m at worst is a good part of what you came to look at.
%
%   psectime_mooring also has a 'sigma' MODE, where it interpolates each
%   Thorpe-sorted profile onto fixed target surfaces so a row follows one
%   ISOPYCNAL instead of one depth, and heave leaves the colour field entirely.
%   That is the better view for separating a water-mass change from a heave, and
%   it is NOT what this function draws -- switch the ztype argument in the four
%   calls below if you want it.
%
%   zlim is passed through as psectime_mooring's zmax: a cut on the sensors'
%   MEAN depth, not an axis limit, so a sensor that gets knocked below zlim for
%   part of the record is still drawn.
%
%   ---------------------------------------------------------------------
%   zlim IS A SCALAR MAX DEPTH, NOT AN AXIS RANGE
%   ---------------------------------------------------------------------
%   The name invites [zmin zmax], because that is what MATLAB's built-in zlim
%   does, and this variable shadows it inside this function. It is not that:
%   psectime_mooring and psecline_mooring both take a single max depth with []
%   meaning no cut. A two-element vector is therefore rejected here rather than
%   passed on, where it would reach a comparison against a matrix and either
%   error somewhere unhelpful or silently keep the wrong rows.
%
% as part of pipeline for computing IT flux on FLEAT moorings
%
% author TEJ, with help from Claude
% 29 September 2026, zlim and figbase 16 October 2026

%% ---- arguments ------------------------------------------------------
if nargin < 2, cfg       = struct(); end
if nargin < 3, temp_cont = []; end
if nargin < 4, sal_cont  = []; end
if nargin < 5, uv_cont   = []; end
if nargin < 6, time_lim  = []; end
if nargin < 7, zlim      = []; end
if nargin < 8 || isempty(figbase), figbase = 0; end

if ~isempty(zlim)
    if ~isnumeric(zlim) || ~isscalar(zlim) || ~isfinite(zlim) || zlim <= 0
        error('psectime_glance_mooring:zlim', ...
            ['zlim is a single max depth in metres, e.g. 300, or [] for the ' ...
             'full column -- not an axis range. Got %s.'], mat2str(zlim));
    end
end

% The black contours drawn over the colour: which field, at what levels, and
% which of those are drawn heavy and labelled.
%
% sig_var IS THE OFF SWITCH, NOT sig_cont. psectime_mooring gates the overlay
% on var2 being non-empty, whereas an empty cont2 or cont2h means AUTOSET, not
% none:
%     cont2  = []  ->  12 levels over the 2nd-98th percentile of the field
%     cont2h = []  ->  every 5th of cont2
% So blanking the two level lists does not remove the contours, it hands the
% choice of levels to the autoset -- which is why sig_cont = [] still drew
% them. Only sig_var = [] gives a clean colour field with nothing over it.
%
% Currently OFF. sig_var = 'sigma' puts the isopycnals back, and then the
% levels below start mattering again: 20:.2:30 and 20:30 were the values used
% up to now, chosen to span every water mass in the array so that one set
% served every mooring.
sig_var   = [];             % 'sigma' to overlay isopycnals, [] for none
sig_cont  = [];             % [] autosets; ignored while sig_var is []
sig2_cont = [];             % the heavy, labelled subset

p = struct('fig1',[], 'fig2',[], 'T',[], 'S',[], 'u',[], 'v',[], 'sigma',[]);

%% ---- figure 1: T, S, u, v against depth -----------------------------
p.fig1 = figure(figbase + 1); clf(p.fig1)
% Four stacked panels need the height, or the colourbars eat the axes.
set(p.fig1, 'Color', 'w', 'Position', [60 60 1100 820])

subplot(4,1,1)
p.T = psectime_mooring(F_h, 'T', sig_var, temp_cont, sig_cont, sig2_cont, 'depth', zlim, cfg, time_lim);
subplot(4,1,2)
p.S = psectime_mooring(F_h, 'S', sig_var, sal_cont,  sig_cont, sig2_cont, 'depth', zlim, cfg, time_lim);
subplot(4,1,3)
p.u = psectime_mooring(F_h, 'u', sig_var, uv_cont,   sig_cont, sig2_cont, 'depth', zlim, cfg, time_lim);
subplot(4,1,4)
p.v = psectime_mooring(F_h, 'v', sig_var, uv_cont,   sig_cont, sig2_cont, 'depth', zlim, cfg, time_lim);

% One zoom drives all four, in time AND in depth. The handles are used rather
% than findall(gcf,'Type','axes') because that would also catch anything else
% living in the figure, and the four panels are exactly what should move
% together.
%
% 'xy' is safe here because all four panels share a vertical coordinate by
% construction: in 'depth' mode psectime_mooring sets Y from the depth matrix,
% which does not depend on which variable is being contoured, and all four get
% the same zlim. So the four YLims are already identical and linking them
% changes nothing about the initial view -- it only ties the zooming. If you
% ever point these calls at different ztypes, or give one panel its own zlim,
% drop back to 'x'.
%
% linkaxes applies the FIRST axes' limits to the rest, so p.T leads. That
% matters only if the panels ever stop agreeing.
linkaxes([p.T.ax p.S.ax p.u.ax p.v.ax], 'xy');

%% ---- figure 2: knockdown, whole record ------------------------------
p.fig2 = figure(figbase + 2); clf(p.fig2)
set(p.fig2, 'Color', 'w')
p.sigma = psecline_mooring(F_h, 'sigma', [], 'depth', zlim, cfg);

end
