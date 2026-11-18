function mode_movie(eta, zref, M, ahat, fit, th, opt)
% MODE_MOVIE  Step through profiles: obs eta vs mode 1, 1+2, 1+2+3, 1+2+3+4.
%
%   mode_movie(eta, zref, M, ahat, fit, th, opt)
%
%   opt .frames    vector of time indices to animate (default: every 3rd hour
%                  over the first 10 days that have a valid fit)
%       .outfile   'modes_F2.mp4' -> writes a movie. Empty ('') just animates
%                  on screen without writing. Default 'modes_F2.mp4'.
%       .fps       frames per second (default 8)
%       .etalim    x limits in metres (default symmetric, 99th pct of |eta|)
%       .pplot     true (default) to add a right-hand panel with p'(z)
%
%   Left panel: observed eta at the instrument levels (black circles) with the
%   cumulative modal reconstructions drawn as smooth full-depth curves. Right
%   panel: the corresponding p' built from the same amplitudes, so you can watch
%   how much of the pressure signal each mode is actually carrying. The title
%   carries the cumulative variance explained, which is your quantitative answer
%   to "how many modes to go to" -- Alford's note leaves that as "know when you
%   look at it", and this turns it into a number.

if nargin < 7, opt = struct(); end
if ~isfield(opt,'outfile'), opt.outfile = 'modes_F2.mp4'; end
if ~isfield(opt,'fps'),     opt.fps     = 8;  end
if ~isfield(opt,'pplot'),   opt.pplot   = true; end

nm = size(ahat,1);
valid = find(all(isfinite(ahat),1));
if isempty(valid), error('mode_movie:nofits','No timesteps have a complete fit.'); end

if ~isfield(opt,'frames') || isempty(opt.frames)
    i0 = valid(1);
    opt.frames = valid(valid >= i0 & valid <= i0 + 240);
    opt.frames = opt.frames(1:3:end);
end
if ~isfield(opt,'etalim') || isempty(opt.etalim)
    L = prctile(abs(eta(:)), 99); opt.etalim = [-L L]*1.15;
end

% full-depth mode shapes and pressure modes on the model grid
Zf = M.Zeta(:,1:nm);
Pf = M.P(:,1:nm);
plim = max(abs(Pf*ahat(:,opt.frames)),[],'all')*1.15;
if ~isfinite(plim) || plim == 0, plim = 1; end

cols = lines(nm);
lbl  = arrayfun(@(m) sprintf('modes 1-%d', m), 1:nm, 'uni', 0);
lbl{1} = 'mode 1';

fh = figure('Position',[100 100 900 720],'Color','w');
% Locked geometry: a figure resized mid-render changes the captured frame
% size, which shears the movie. See capture_frame.
set(fh, 'Resize', 'off', 'PaperPositionMode', 'auto');
FRAME_DPI = 150;   % capture resolution, dots per inch
fsz = [];          % frame size, locked by the first capture_frame call
if ~isempty(opt.outfile)
    vw = VideoWriter(opt.outfile,'MPEG-4'); vw.FrameRate = opt.fps; open(vw);
end

for k = 1:numel(opt.frames)
    it = opt.frames(k);
    clf(fh);

    % ---------------- left: eta ----------------
    if opt.pplot, subplot(1,2,1); end
    hold on; box on
    h = gobjects(nm+1,1);
    for m = 1:nm
        h(m) = plot(Zf(:,1:m)*ahat(1:m,it), M.z, '-', 'Color', cols(m,:), 'LineWidth', 1.6);
    end
    h(nm+1) = plot(eta(:,it), zref, 'ko', 'MarkerSize', 5, ...
                   'MarkerFaceColor',[.85 .85 .85], 'LineWidth', 0.8);
    xline(0,'k:');
    set(gca,'YLim',[0 M.H],'XLim',opt.etalim);
    xlabel('\eta (m)'); ylabel('z (m)');
    legend(h, [lbl {'observed'}], 'Location','southeast','Box','off','FontSize',8);
    set(gca, 'YDir', 'reverse');
    ve = fit.varexp(:,it);
    title(sprintf('%s\nvar. explained  %s', datestr(th(it),'dd mmm yyyy HH:MM'), ...
          strjoin(arrayfun(@(v) sprintf('%.2f',v), ve(:).','uni',0), ' / ')), ...
          'FontWeight','normal','FontSize',9);

    % ---------------- right: p' ----------------
    if opt.pplot
        subplot(1,2,2); hold on; box on
        for m = 1:nm
            plot(Pf(:,1:m)*ahat(1:m,it), M.z, '-', 'Color', cols(m,:), 'LineWidth', 1.6);
        end
        xline(0,'k:');
        set(gca,'YLim',[0 M.H],'XLim',[-plim plim]);
        xlabel('p'' [Pa]'); ylabel('z [m]');
        title(sprintf('amplitudes [m]: %s', ...
              strjoin(arrayfun(@(a) sprintf('%+.1f',a), ahat(:,it).','uni',0), '  ')), ...
              'FontWeight','normal','FontSize',9);
    end
    set(gca, 'YDir', 'reverse');
    drawnow;
    if ~isempty(opt.outfile)
        [im, fsz] = capture_frame(fh, FRAME_DPI, fsz);
        writeVideo(vw, im);
    end
end

if ~isempty(opt.outfile)
    close(vw);
    fprintf('wrote %s  (%d frames)\n', opt.outfile, numel(opt.frames));
end
end