function pmode_sweep(S, maxc, ttl, fignum)
%PMODE_SWEEP  The "how many modes do these levels support?" figure.
%
%   pmode_sweep(S, maxc, ttl, fignum)
%
%   S  .m          [1 x nmax]   number of modes fit
%      .pct        [nmax x 1]   cumulative variance explained, NO gate
%      .pctc       [nmax x 1]   the same, scored only on hours passing maxc
%      .frac       [nmax x 1]   fraction of hours passing maxc
%      .cond       [nmax x 3]   cond(A_w): p10, median, p90
%      .pct_extra  [nmax x 1]   optional second ungated curve, drawn grey
%      .pct_name / .pctc_name / .extra_name   legend text
%   maxc    the conditioning limit the gate used. Inf is legal (gate off).
%   ttl     figure title
%   fignum  figure number
%
%   WHY TWO RULERS
%   --------------
%   The quantities do not share a range and cannot share an axis. Variance
%   explained has a CEILING of 1 and no floor -- a reconstruction can be
%   arbitrarily worse than predicting zero, and F1 printed -1.9 at m = 1.
%   cond has a FLOOR of 1 (perfectly distinguishable columns) and no ceiling;
%   F1 reaches 1e4 by m = 8. So cond goes on the right, log-scaled.
%
%   HOW TO READ IT
%   --------------
%   Both of the interesting curves RISE with m, which is the whole point:
%   variance explained rises because more modes always fit better, and cond
%   rises because more modes are always harder to tell apart. Dropping a
%   column cannot raise the largest singular value nor lower the smallest, so
%   cond(A(:,1:m)) <= cond(A(:,1:m+1)) always -- the cost curve is monotone by
%   construction, not by luck.
%
%   Read the crossing: the largest m whose cond still clears maxcond on most
%   hours. The gap between the gated and ungated variance curves IS the
%   conditioning. The .frac curve is there because a gated curve computed on
%   3% of the hours is not a result, and nothing else on the figure says so.
%
%   author TEJ, with help from Claude
%   10 November 2026

if nargin < 4 || isempty(fignum), fignum = 9; end
if nargin < 3, ttl = ''; end
if ~isfield(S,'pct_name'),   S.pct_name   = 'varexp, no gate';     end
if ~isfield(S,'pctc_name'),  S.pctc_name  = 'varexp, gated';       end
if ~isfield(S,'extra_name'), S.extra_name = 'raw varexp, no gate'; end

m    = S.m(:)';
nmax = numel(m);

% Colours are explicit: yyaxis gives each side its own ruler colour and will
% not hand four left-hand lines four distinguishable ones.
cNG = [0.00 0.45 0.74];   cG = [0.85 0.33 0.10];
cFR = [0.47 0.67 0.19];   cC = [0.49 0.18 0.56];

figure(fignum); clf; set(gcf, 'Color', 'w')

yyaxis left
  h1 = plot(m, S.pct(:),  'o-',  'LineWidth',1.8, 'Color',cNG); hold on
  h2 = plot(m, S.pctc(:), 's--', 'LineWidth',1.8, 'Color',cG);
  h3 = plot(m, S.frac(:), '^:',  'LineWidth',1.2, 'Color',cFR);
  hh  = [h1 h2 h3];
  lab = {S.pct_name, S.pctc_name, 'hours passing the gate'};
  if isfield(S,'pct_extra') && ~isempty(S.pct_extra)
      h4 = plot(m, S.pct_extra(:), '-', 'LineWidth',0.9, 'Color',[.55 .55 .55]);
      hh(end+1) = h4;  lab{end+1} = S.extra_name;
  end
  yline(0, 'k-', 'LineWidth',0.5);
  ylabel('fraction of variance  /  fraction of hours')
  lo = min([S.pct(:); S.pctc(:); 0]);
  if ~isfinite(lo), lo = 0; end             % sweep starved: everything NaN
  ylim([lo - 0.1, 1.05])
  set(gca, 'YColor', 'k')

yyaxis right
  h5 = plot(m, S.cond(:,2), 'o-', 'LineWidth',1.8, 'Color',cC); hold on
  plot(m, S.cond(:,1), ':', 'LineWidth',1.0, 'Color',cC)
  plot(m, S.cond(:,3), ':', 'LineWidth',1.0, 'Color',cC)
  % Inf is a legal maxcond (the gate off, for a diagnostic run) and yline
  % rejects a non-finite value. With the gate off the gated and ungated curves
  % coincide, which is the honest picture anyway.
  if isfinite(maxc)
      yline(maxc, '--', sprintf('maxcond = %g', maxc), 'LineWidth',1.4);
  end
  if any(isfinite(S.cond(:))), set(gca, 'YScale', 'log'); end
  ylabel('cond(A_w)   [log, 1 is ideal]')
  set(gca, 'YColor', cC)

xlim([0.5 nmax+0.5]); xticks(m); grid on
xlabel('number of modes fit, m')
legend([hh h5], [lab {'cond(A_w), median (10-90% dotted)'}], ...
       'Location','southoutside', 'Orientation','vertical', 'Box','off')
title(ttl)
end
