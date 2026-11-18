function S = uv_mode_sweep(Mh, Mx, opt)
%UV_MODE_SWEEP  How many VELOCITY modes do the ADCP bins support?
%
%   S = uv_mode_sweep(Mh, Mx, opt)
%
%   Mh         the hourly mooring struct: needs .u, .v, .depth
%   Mx         a modes struct from build_modes with AT LEAST opt.nmax modes,
%              so Mx.Zu0 has nmax+1 columns. The nbc-mode struct the pipeline
%              fits with is not enough when nmax > nbc -- build a wider one.
%   opt .nmax     how many baroclinic modes to sweep to
%       .w        [ns x 1] the weights fit_uv_modes actually used (U.w)
%       .maxcond  the gate, matching cfg.fit_maxcond
%       .minobs   minimum good bins in an hour
%
%   S  .m .pct .pctc .frac .cond .nhr   -- the inputs pmode_sweep wants.
%
%   WHY THE VELOCITY SIDE NEEDS ITS OWN SWEEP
%   -----------------------------------------
%   It is a different matrix, and at F1 it is the BINDING one. Three
%   differences from the displacement fit, all of which make it harder:
%
%     1. It solves for nm+1 columns, not nm. Column 1 is the barotropic mode 0
%        (fit_uv_modes.m:156), so asking for 4 baroclinic modes is a 5-column
%        problem.
%     2. The ADCP spans less of the column than the thermistor string.
%     3. The design matrix is REBUILT EVERY TIMESTEP, because the sensors move
%        -- 56 m of median knockdown at F2. So cond varies hour to hour from
%        both which bins survive AND where they were, whereas the displacement
%        fit works on a fixed zref and varies only through dropout. Expect a
%        wider 10-90% band here than on the eta figure.
%
%   The hour gate is max(minobs, nmax+3) for EVERY m, so all m are scored on
%   the same hours and the curve is a comparison rather than a mix of
%   subsamples. Same convention as the eta sweep in fleat_IT_flux cell 6.
%
%   WHAT FIELD THIS SCORES, AND WHAT IT IS COMPARABLE TO
%   ----------------------------------------------------
%   The UNFILTERED velocity, with only the time-mean profile removed -- which
%   is what fit_uv_modes itself fits, because the pipeline fits first and
%   bandpasses the amplitudes afterwards. There is no band-passed u,v at the
%   bins to score against.
%
%   So S.pct is the analogue of the displacement sweep's RAW eta curve (the
%   grey one), NOT of its band curve. Do not read the two figures' main curves
%   against each other: at F1 this lands at 0.54 for m = 1, which sits with
%   raw eta's 0.51 and nowhere near band eta's 0.70. Like the raw curve, it
%   includes the broadband continuum that m modes were never going to
%   represent, so it is a floor on what the modes do for the tide, not a
%   measure of it.
%
%   author TEJ, with help from Claude
%   10 November 2026

if ~isfield(opt,'maxcond') || isempty(opt.maxcond), opt.maxcond = 100; end
if ~isfield(opt,'minobs')  || isempty(opt.minobs),  opt.minobs  = 8;   end

nmax = opt.nmax;
w    = opt.w(:);
maxc = opt.maxcond;

if ~isfield(Mx,'Zu0')
    error('uv_mode_sweep:stalemodes', 'Mx has no Zu0 field -- re-run build_modes.');
end
if size(Mx.Zu0,2) < nmax+1
    error('uv_mode_sweep:nmodes', ...
        ['Mx.Zu0 has %d columns, need %d (mode 0 plus %d baroclinic). Build ' ...
         'the sweep modes with nmax, not with cfg.nbc.'], ...
        size(Mx.Zu0,2), nmax+1, nmax);
end

u = Mh.u;  v = Mh.v;  Dz = Mh.depth;
[ns, nt] = size(u);
if numel(w) ~= ns
    error('uv_mode_sweep:w', 'opt.w has %d entries but u has %d rows.', numel(w), ns);
end

% Same first step as the fit: remove the time-mean profile, so the sweep scores
% the same field fit_uv_modes scores.
up0 = u - mean(u, 2, 'omitnan');
vp0 = v - mean(v, 2, 'omitnan');

Zall = Mx.Zu0(:,1:nmax+1);
zlo  = Mx.z(1);  zhi = Mx.z(end);
need = max(opt.minobs, nmax+3);

cn   = nan(nt, nmax);
num  = zeros(nmax,1);  den  = 0;
numg = zeros(nmax,1);  deng = zeros(nmax,1);
nhr  = 0;

for it = 1:nt
    zi = Dz(:,it);
    ok = isfinite(zi) & isfinite(up0(:,it)) & isfinite(vp0(:,it));
    if nnz(ok) < need, continue; end

    A = interp1(Mx.z, Zall, min(max(zi(ok), zlo), zhi));
    if any(~isfinite(A(:))), continue; end

    Wk = sqrt(w(ok));
    Aw = A .* Wk;
    eu = up0(ok,it);  ev = vp0(ok,it);
    v0 = sum(w(ok).*(eu.^2 + ev.^2));
    if ~(v0 > 0), continue; end

    nhr = nhr + 1;
    den = den + v0;
    for m = 1:nmax
        Am = Aw(:,1:m+1);                   % +1: mode 0 is always in the fit
        cn(it,m) = cond(Am);
        au = Am \ (eu.*Wk);
        av = Am \ (ev.*Wk);
        ru = eu - A(:,1:m+1)*au;
        rv = ev - A(:,1:m+1)*av;
        r  = sum(w(ok).*(ru.^2 + rv.^2));   % u and v scored together
        num(m) = num(m) + r;
        if cn(it,m) <= maxc
            numg(m) = numg(m) + r;  deng(m) = deng(m) + v0;
        end
    end
end

S.m    = 1:nmax;
S.nhr  = nhr;
S.pct  = nan(nmax,1);
S.pctc = nan(nmax,1);
S.cond = nan(nmax,3);
S.frac = nan(nmax,1);

if den > 0
    S.pct = 1 - num/den;
end
g = deng > 0;                               % NaN where no hour passed the gate
S.pctc(g) = 1 - numg(g)./deng(g);

scored = any(isfinite(cn), 2);
if any(scored)
    C = cn(scored,:);
    S.cond = [prctile(C,10,1); median(C,1); prctile(C,90,1)]';
    S.frac = mean(C <= maxc, 1)';
end
end
