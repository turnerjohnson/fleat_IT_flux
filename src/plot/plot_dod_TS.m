function out = plot_dod_TS(opt)
% PLOT_DOD_TS  T-S diagram of the real M1-M5 instrument data, with isopycnals.
%
%   out = plot_dod_TS()
%   out = plot_dod_TS(opt)
%
%   opt .datadir   where the M*_TSPbar.mat live
%       .lat .lon  for the TEOS-10 conversions. Default 7.9 N / 134.6 E --
%                  PLACEHOLDERS, same caveat as build_dod_mooring
%       .stride    subsample in time. Default 6, i.e. hourly from 10-minute
%       .axes      'CT-SA' (default) | 't-SP'   which pair to plot; see below
%       .qc        true (default) show the rejected points separately
%
%   out .SA .CT .p .m   the plotted points, so a driver can re-use them
%       .stats          per-mooring, per-sensor summary (the table view)
%       .offset         how far the TEOS-10 conversion moved the data
%
%   WHICH TEMPERATURE AND SALINITY THESE ACTUALLY ARE
%   The files hold IN-SITU temperature and PRACTICAL salinity -- TSdbar in degC
%   and Sdbar in psu. build_dod_mooring stores exactly that in F.T and F.S; it
%   converts to SA and CT only inside the N2 calculation and does not keep
%   them. So no, nothing upstream of here is in conservative temperature.
%
%   This function converts properly before plotting, because a T-S diagram with
%   isopycnals has to be internally consistent:
%       SA  = gsw_SA_from_SP(SP, p, lon, lat)      g/kg, absolute salinity
%       CT  = gsw_CT_from_t(SA, t, p)              degC, conservative temp
%       sig0 = gsw_sigma0(SA, CT)                  kg/m3, ref 0 dbar
%   Pressure comes straight from PSdbar, which already IS pressure in dbar, so
%   nothing is reconstructed from a nominal depth.
%
%   opt.axes = 't-SP' plots the raw in-situ pair instead, for comparison. The
%   isopycnals are still computed from SA and CT -- there is no correct way to
%   contour potential density on in-situ axes, so that mode labels them as
%   indicative and the printed offsets say how far the axes have shifted.
%
%   ONE ASSUMPTION THAT IS NOT CHECKABLE FROM THE FILES
%   Nothing in these .mat files records the temperature scale or the salinity
%   definition. In-situ ITS-90 and practical salinity is the overwhelmingly
%   likely reading for a 2016 mooring product, and the values are consistent
%   with it, but it is an inference. If the product turns out to be CT and SA
%   already then this function double-converts, which would show up as a ~0.1
%   degC and ~0.17 g/kg shift -- the printed offsets are the thing to compare
%   against if anyone ever finds the processing notes.
%
%   WHY THE MOORINGS ARE FACETED RATHER THAN COLOUR-CODED
%   A scatter puts every pair of series on screen at once, and the validated
%   categorical palette only clears the all-pairs colourblind floors for its
%   first three slots. Five moorings in five hues would fail that gate, so
%   mooring identity is carried by POSITION -- one small panel each -- and
%   colour is freed for pressure, which is a magnitude and belongs on a
%   sequential ramp anyway. The three discrete marks that do carry identity
%   (isopycnals, the S(T) curve, the rejected points) were validated all-pairs:
%   worst normal-vision dE 31.6, worst CVD dE 23.8, against floors of 15 and 8.
%
%   WHAT TO LOOK FOR
%   The top panel is every mooring pooled. The tight lower-left hook is the
%   deep water -- that is the part where S varies by 0.004-0.011 psu and why
%   build_dod_mooring is willing to hold S constant below the deepest sensor.
%   The spread at the warm end is the near-surface variability that the 'ts'
%   S-method cannot represent, since there one T maps to many S. The black
%   curve is exactly the S(T) relation that opt.s_method 'ts' would impose: if
%   it runs through the middle of a wide cloud rather than along a narrow one,
%   that is the method's error, drawn to scale.

%% ---- arguments ------------------------------------------------------
if nargin < 1 || ~isstruct(opt), opt = struct(); end
if ~isfield(opt,'datadir') || isempty(opt.datadir)
    opt.datadir = '/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/';
end
if ~isfield(opt,'lat'),    opt.lat    = 7.9;      end
if ~isfield(opt,'lon'),    opt.lon    = 134.6;    end
if ~isfield(opt,'stride'), opt.stride = 6;        end
if ~isfield(opt,'axes'),   opt.axes   = 'CT-SA';  end
if ~isfield(opt,'qc'),     opt.qc     = true;     end

if isempty(which('gsw_sigma0'))
    error('plot_dod_TS:nogsw', ...
        ['Needs TEOS-10. addpath(''/Users/tej/Documents/MATLAB/gsw''). ' ...
         'There is no point drawing isopycnals from an approximation here.']);
end

% palette -- validated; see WHY THE MOORINGS ARE FACETED above
C.surface  = [252 252 251]/255;
C.ink      = [ 11  11  11]/255;    % the S(T) curve, and primary text
C.ink2     = [ 82  81  78]/255;    % isopycnals, axes, secondary text
C.critical = [208  59  59]/255;    % qc-rejected points (always labelled)
C.ramp     = [205 226 251; 158 197 244; 109 167 236; 57 135 229; ...
               37 106 191;  24  92 171; 13  54 107]/255;   % blue 100->700
cmap = interp1(linspace(0,1,size(C.ramp,1)), C.ramp, linspace(0,1,256));

%% ---- load the real instrument data ---------------------------------
% The gridded TZdepth product has no measured salinity, so the only real T-S
% pairs in this dataset are the instrument records in TSPbar.
lim = struct('S', [32 35.6], 'T', [3 32], 'P', [10 1200]);
[SP, t, p, mm, ss] = deal([]);
[SPx, tx, px, mmx]  = deal([]);                 % the qc-rejected ones
for m = 1:5
    fb = fullfile(opt.datadir, sprintf('M%d_TSPbar.mat', m));
    if ~isfile(fb)
        warning('plot_dod_TS:missing', 'No file for M%d, skipped.', m);  continue
    end
    B  = load(fb, 'Sdbar', 'TSdbar', 'PSdbar');
    k  = 1:opt.stride:size(B.Sdbar, 1);
    S0 = double(B.Sdbar(k,:));  T0 = double(B.TSdbar(k,:));  P0 = double(B.PSdbar(k,:));
    ns = size(S0, 2);
    sensor = repmat(1:ns, numel(k), 1);
    moor   = m*ones(numel(k), ns);

    good = isfinite(S0) & isfinite(T0) & isfinite(P0) ...
         & S0 >= lim.S(1) & S0 <= lim.S(2) ...
         & T0 >= lim.T(1) & T0 <= lim.T(2) ...
         & P0 >= lim.P(1) & P0 <= lim.P(2);
    bad  = isfinite(S0) & isfinite(T0) & isfinite(P0) & ~good;

    SP = [SP; S0(good)];  t = [t; T0(good)];  p = [p; P0(good)];      %#ok<AGROW>
    mm = [mm; moor(good)]; ss = [ss; sensor(good)];                   %#ok<AGROW>
    SPx = [SPx; S0(bad)]; tx = [tx; T0(bad)]; px = [px; P0(bad)];     %#ok<AGROW>
    mmx = [mmx; moor(bad)];                                           %#ok<AGROW>
end
if isempty(SP), error('plot_dod_TS:nodata', 'No usable T/S pairs found.'); end

%% ---- TEOS-10 ------------------------------------------------------
SA  = gsw_SA_from_SP(SP, p, opt.lon, opt.lat);
CT  = gsw_CT_from_t(SA, t, p);
SAx = gsw_SA_from_SP(SPx, px, opt.lon, opt.lat);
CTx = gsw_CT_from_t(SAx, tx, px);

out.offset.SA_minus_SP = mean(SA - SP);
out.offset.CT_minus_t  = mean(CT - t);
fprintf(['plot_dod_TS: %d points (every %dth sample), %d qc-rejected\n' ...
         '  TEOS-10 shift: SA - SP = %+.4f g/kg mean, CT - t = %+.4f degC mean\n'], ...
    numel(SA), opt.stride, numel(SAx), out.offset.SA_minus_SP, out.offset.CT_minus_t);

usingCT = strcmpi(opt.axes, 'CT-SA');
if usingCT
    X = SA;  Y = CT;  Xx = SAx;  Yx = CTx;
    xl = 'absolute salinity S_A  (g kg^{-1})';  yl = 'conservative temperature \Theta  (\circC)';
else
    X = SP;  Y = t;   Xx = SPx;  Yx = tx;
    xl = 'practical salinity S_P  (psu)';       yl = 'in-situ temperature t  (\circC)';
end

%% ---- isopycnal grid ------------------------------------------------
% Always computed from SA and CT, whichever axes are drawn.
padx = 0.04*range(X);  pady = 0.04*range(Y);
xv = linspace(min(X)-padx, max(X)+padx, 160);
yv = linspace(min(Y)-pady, max(Y)+pady, 160);
[XG, YG] = meshgrid(xv, yv);
if usingCT
    SIG = gsw_sigma0(XG, YG);
else
    % map the in-situ axes to SA/CT at a representative pressure so the
    % contours are at least self-consistent; they are indicative only
    pref = median(p);
    SAg  = gsw_SA_from_SP(XG, pref, opt.lon, opt.lat);
    SIG  = gsw_sigma0(SAg, gsw_CT_from_t(SAg, YG, pref));
end

%% ---- figure --------------------------------------------------------
f = figure('Color', C.surface, 'Position', [60 40 1120 980]);
tl = tiledlayout(f, 3, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

% -- pooled panel ----------------------------------------------------
ax = nexttile(tl, 1, [2 5]);
hold(ax, 'on');  set(ax, 'Color', C.surface, 'XColor', C.ink2, 'YColor', C.ink2, ...
    'FontSize', 11, 'Layer', 'top', 'TickDir', 'out', 'Box', 'off');

[CS, hC] = contour(ax, xv, yv, SIG, 14, 'Color', C.ink2, 'LineWidth', 0.5);
clabel(CS, hC, 'FontSize', 8, 'Color', C.ink2, 'LabelSpacing', 400);

sc = scatter(ax, X, Y, 4, p, 'filled');
sc.MarkerFaceAlpha = 0.06;                      % a density cloud, not marks
colormap(ax, cmap);  set(ax, 'CLim', [0 max(p)]);
cb = colorbar(ax);  cb.Direction = 'reverse';   % deep at the bottom
cb.Label.String = 'pressure (dbar)';  cb.Color = C.ink2;  cb.Label.Color = C.ink;

% the S(T) curve that opt.s_method 'ts' would impose, built the same way
[tc, sc_, rel] = ts_relation(t, SP);
if usingCT
    pc  = median(p)*ones(size(tc));
    SAc = gsw_SA_from_SP(sc_, pc, opt.lon, opt.lat);
    Xc  = SAc;  Yc = gsw_CT_from_t(SAc, tc, pc);
else
    Xc  = sc_;  Yc = tc;
end

% Draw the curve in two styles. Solid where the bins behind it hold real
% water; dashed across the void between the 4th and deepest instruments,
% where the curve is a straight line between two unrelated water masses and
% not a measured relation. Smoothing the void away would have hidden it --
% ts_relation marks it instead, and s_method 'ts' leaves S NaN there.
hFit = plot(ax, Xc, Yc, '-', 'Color', C.ink, 'LineWidth', 2);
leg = hFit;  lbl = {'S(T) relation used by s\_method ''ts'''};
for r = 1:size(rel.gap_range, 1)
    lo = rel.gap_range(r,1);  hi = rel.gap_range(r,2);
    kk = find(tc <= lo, 1, 'last');  k2 = find(tc >= hi, 1, 'first');
    if ~isempty(kk) && ~isempty(k2)
        plot(ax, Xc([kk k2]), Yc([kk k2]), '-', 'Color', C.surface, 'LineWidth', 3.5);
        hGap = plot(ax, Xc([kk k2]), Yc([kk k2]), ':', 'Color', C.critical, ...
            'LineWidth', 2.5);
        if r == 1
            leg(end+1) = hGap;  %#ok<AGROW>
            lbl{end+1} = sprintf(['no data %.1f-%.1f \\circC: unsupported, ' ...
                'S left NaN'], lo, hi);  %#ok<AGROW>
        end
    end
end
if opt.qc && ~isempty(Xx)
    hBad = plot(ax, Xx, Yx, 'x', 'Color', C.critical, 'MarkerSize', 7, 'LineWidth', 1.2);
    leg(end+1) = hBad;
    lbl{end+1} = sprintf('qc-rejected (%d points, S outside [%g %g])', ...
        numel(Xx), lim.S(1), lim.S(2));
end
legend(ax, leg, lbl, 'Location', 'southeast', 'TextColor', C.ink, ...
    'Color', C.surface, 'EdgeColor', C.ink2, 'FontSize', 10);
xlabel(ax, xl, 'Color', C.ink);  ylabel(ax, yl, 'Color', C.ink);
title(ax, 'FLEAT DOD moorings M1-M5: measured T-S, all sensors pooled', ...
    'Color', C.ink, 'FontSize', 13, 'FontWeight', 'bold');
subtitle(ax, sprintf(['contours \\sigma_0 (kg m^{-3}); %d points, every %dth sample; ' ...
    'isopycnals from S_A and \\Theta'], numel(X), opt.stride), ...
    'Color', C.ink2, 'FontSize', 10);
grid(ax, 'off');

% -- small multiples, one per mooring --------------------------------
for m = 1:5
    axm = nexttile(tl, 10 + m);
    hold(axm, 'on');  set(axm, 'Color', C.surface, 'XColor', C.ink2, ...
        'YColor', C.ink2, 'FontSize', 9, 'TickDir', 'out', 'Box', 'off');
    contour(axm, xv, yv, SIG, 14, 'Color', C.ink2, 'LineWidth', 0.3);
    g = mm == m;
    s2 = scatter(axm, X(g), Y(g), 3, p(g), 'filled');
    s2.MarkerFaceAlpha = 0.10;
    colormap(axm, cmap);  set(axm, 'CLim', [0 max(p)]);
    xlim(axm, [xv(1) xv(end)]);  ylim(axm, [yv(1) yv(end)]);
    title(axm, sprintf('M%d  (%d sensors)', m, numel(unique(ss(g)))), ...
        'Color', C.ink, 'FontSize', 10);
    if m == 1, ylabel(axm, yl, 'Color', C.ink, 'FontSize', 9); end
    if m == 3, xlabel(axm, xl, 'Color', C.ink, 'FontSize', 9); end
    if m > 1,  set(axm, 'YTickLabel', []); end
end

out.SA = SA;  out.CT = CT;  out.p = p;  out.m = mm;  out.sensor = ss;
out.ts_curve = struct('t', tc, 'S', sc_, 'info', rel);
out.fig = f;

fprintf(['  S(T) curve: %d equal-count bins, %d-%d points each, ' ...
         '%d void(s) flagged\n'], rel.nbin_used, min(rel.n), max(rel.n), ...
         size(rel.gap_range, 1));
for r = 1:size(rel.gap_range, 1)
    fprintf('    void %.2f - %.2f degC (no instrument samples this water)\n', ...
        rel.gap_range(r,1), rel.gap_range(r,2));
end

%% ---- the table view -----------------------------------------------
% Required alongside the figure: the same content as numbers, so nothing here
% depends on reading a colour off a dense cloud.
fprintf('\n  per-sensor summary (the table view)\n');
fprintf('  moor sens   p med    SA med    CT med    SA std    CT std   sig0 med\n');
k = 0;
for m = 1:5
    for s = unique(ss(mm == m))'
        g = mm == m & ss == s;  k = k + 1;
        out.stats(k,:) = [m s median(p(g)) median(SA(g)) median(CT(g)) ...
                          std(SA(g)) std(CT(g)) median(gsw_sigma0(SA(g), CT(g)))];
        fprintf('   M%d   %d  %7.1f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n', ...
            out.stats(k,:));
    end
end
end
