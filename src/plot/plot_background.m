function h = plot_background(B, cfg, fignum)
%PLOT_BACKGROUND  The three-panel look at the mean profile, N^2 and T_b.
%
%   h = plot_background(B, cfg)
%   h = plot_background(B, cfg, fignum)
%
%   B        from f_background_v2
%   cfg      mooring config; only cfg.mooring and cfg.mld are read
%   fignum   figure number to reuse and clear; omit or [] for a new figure
%
%   Panel 1  mean sigma_0, full grid and the binned reference levels
%   Panel 2  BOTH N^2 estimates. The one currently selected is drawn heavy and
%            named in the legend; the other is drawn light so you can see what
%            you would be switching to. Grey dots are the raw gridded median
%            before smoothing.
%   Panel 3  buoyancy period of the selected profile, as a sanity check
%
%   This was the showPlots block sitting inline in both drivers, where F2 and
%   F4 had two copies drifting apart. Read the middle panel, then switch with
%   the commented SET_N2_SOURCE line just below the call if you want the other
%   one.

mname = 'mooring';
if isfield(cfg,'mooring') && ~isempty(cfg.mooring), mname = cfg.mooring; end
sel   = 'sigma';
if isfield(B,'N2_source'), sel = B.N2_source; end
hasg  = isfield(B,'N2full_grid') && ~isempty(B.N2full_grid);

if nargin < 3 || isempty(fignum)
    h.fig = figure('Color','w');
else
    h.fig = figure(fignum); clf(h.fig); set(h.fig,'Color','w')
end

%% ---- 1. mean density ------------------------------------------------
h.ax(1) = subplot(1,3,1);
plot(B.sig_full, B.zfull, 'LineWidth',1.3); hold on
plot(B.sig_t, B.zref, 'o', 'MarkerSize',4)
if isfield(cfg,'mld') && ~isempty(cfg.mld), yline(cfg.mld,'r--','MLD'); end
xlabel('\sigma_\theta'); ylabel('z [m]')
title('mean profile + mixed layer')
legend('full','binned','Location','southwest','Box','off')
set(gca,'YDir','reverse')

%% ---- 2. the two N^2 estimates --------------------------------------
h.ax(2) = subplot(1,3,2);
lw = @(name) 2.0*strcmp(sel,name) + 0.8*~strcmp(sel,name);
lab = {};
semilogx(B.N2full_sigma, B.zfull, 'LineWidth', lw('sigma')); hold on
lab{end+1} = ptr('from \sigma', sel, 'sigma');
if hasg
    semilogx(B.N2full_grid, B.zfull, 'LineWidth', lw('gridded'))
    lab{end+1} = ptr('F.N2 median', sel, 'gridded');
    semilogx(B.N2grid, B.zgrid, '.', 'MarkerSize',8, 'Color',[.55 .55 .55])
    lab{end+1} = 'F.N2 median, raw';
    yline(B.zgrid(1), ':', 'top of F.N2', 'LabelHorizontalAlignment','left');
end
xlabel('N^2 [s^{-2}]'); title(sprintf('N^2  (%s: %s)', mname, sel))
legend(lab, 'Location','southeast', 'Box','off')
set(gca,'YDir','reverse')

%% ---- 3. sanity check ------------------------------------------------
h.ax(3) = subplot(1,3,3);
plot(2*pi./sqrt(B.N2full)/60, B.zfull, 'LineWidth',1.3)
xlabel('buoyancy period [min]'); title('sanity check')
set(gca,'YDir','reverse')

linkaxes(h.ax, 'y')
end

% =========================================================================
function s = ptr(base, sel, name)
%PTR  Mark the selected curve in the legend so the figure says what is in use.
if strcmp(sel, name)
    s = [base '  <-- in use'];
else
    s = base;
end
end
