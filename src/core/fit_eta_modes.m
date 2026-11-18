function [ahat, fit] = fit_eta_modes(eta, zref, M, opt)
% FIT_ETA_MODES  Least-squares fit of displacement modes to moored eta, timestep by timestep.
%
%   [ahat, fit] = fit_eta_modes(eta, zref, M, opt)
%
%   eta   [nz x nt]  isopycnal displacement, POSITIVE UP, NaNs allowed.
%   zref  [nz x 1]   depth label of each eta row (positive down).
%   M                struct from build_modes.
%   opt   .nmodes    how many modes to solve for (default: all in M.Zeta)
%         .weight    'wkb' (default) | 'none'
%         .minobs    minimum good levels required to attempt a fit (default 8)
%         .maxcond   reject a timestep whose weighted design matrix is worse
%                    conditioned than this (default 100, Inf disables)
%
%   ahat  [nmodes x nt]  modal displacement amplitudes in metres.
%   fit   .recon   [nz x nt x nmodes]  cumulative reconstruction at zref:
%                                       (:,:,m) uses modes 1..m
%         .resid   [nz x nt]  eta minus the full nmodes reconstruction
%         .varexp  [nmodes x nt]  fraction of eta variance explained by modes 1..m
%         .ngood   [1 x nt]   number of good levels used at each timestep
%         .cond    [1 x nt]   condition number of the weighted design matrix
%         .w       [nz x 1]   the weights actually used
%
%   WHY WEIGHTING MATTERS HERE
%   --------------------------
%   Nash et al. (2005) sec. 5 solves a WEIGHTED least-squares problem and defines
%   an "ideally situated" mooring as one with instruments evenly spaced in the
%   WKB-stretched coordinate z' = int_0^z N(z*)/Nbar dz*. F2 is a long way from
%   that. Your bottom eleven RBRs (1389-1503 m, spaced 10 m) sit in weakly
%   stratified water and together occupy about 2 percent of the stretched water
%   column, while single sensors in the thermocline each represent far more of
%   it. An unweighted fit lets that dense bottom cluster pull the solution.
%
%   Each level is weighted by its share of the stretched column. Run with
%   opt.weight = 'none' once to see how much this changes the amplitudes; if the
%   answer is "a lot", the weighting is doing necessary work and should be
%   mentioned when you present.
%
%   WHY minobs IS NOT ENOUGH
%   ------------------------
%   minobs counts levels. It does not ask whether they SPAN the column, and a
%   count test passes happily on levels that cannot separate the modes. At F4
%   there are hours when the whole thermocline drops out and 19 levels survive,
%   all of them below 200 m. Over that range modes 1-4 are nearly collinear:
%   the design matrix reaches cond ~ 6e3 and least squares answers with
%   amplitudes of 20000 to 39000 m that cancel to reproduce a 30 m signal. The
%   reconstruction looks fine and the residual is small, so nothing downstream
%   notices -- but one such hour in a thousand dominates std(ahat) and every
%   amplitude plot drawn from it.
%
%   The conditioning test is the honest version of the count test: it asks
%   whether THESE levels can constrain THESE modes, and returns NaN when they
%   cannot. NaN is the right answer; a 39000 m amplitude is not.

if nargin < 4, opt = struct(); end
if ~isfield(opt,'nmodes') || isempty(opt.nmodes), opt.nmodes = size(M.Zeta,2); end
if ~isfield(opt,'weight') || isempty(opt.weight), opt.weight = 'wkb'; end
if ~isfield(opt,'minobs') || isempty(opt.minobs), opt.minobs = 8;    end
if ~isfield(opt,'maxcond')|| isempty(opt.maxcond),opt.maxcond = 100;  end

zref = zref(:);
[nz, nt] = size(eta);
nm = opt.nmodes;
if numel(zref) ~= nz
    error('fit_eta_modes:size','zref has %d entries but eta has %d rows.', numel(zref), nz);
end
if nm > size(M.Zeta,2)
    error('fit_eta_modes:nmodes','Asked for %d modes but M.Zeta only has %d.', nm, size(M.Zeta,2));
end

% --- design matrix: modes evaluated at the eta reference depths -------------
A = interp1(M.z, M.Zeta(:,1:nm), zref);
if any(~isfinite(A(:)))
    error('fit_eta_modes:interp', ...
        'Some zref fall outside the mode grid (0 to %.0f m). Trim zref.', M.H);
end

% --- weights ---------------------------------------------------------------
switch lower(opt.weight)
    case 'none'
        w = ones(nz,1);
    case 'wkb'
        % each level is weighted by the share of the WKB-stretched water column
        % it represents, i.e. by its spacing in z'
        Nz   = sqrt(M.N2);
        Nbar = trapz(M.z, Nz)/M.H;
        zp   = cumtrapz(M.z, Nz/Nbar);
        zpr  = interp1(M.z, zp, zref);
        edg  = [zpr(1); (zpr(1:end-1)+zpr(2:end))/2; zpr(end)];
        w    = diff(edg);
        w    = w/mean(w);
    otherwise
        error('fit_eta_modes:weight','Unknown weight option "%s".', opt.weight);
end

% --- the fit ---------------------------------------------------------------
ahat  = nan(nm, nt);
ngood = zeros(1, nt);
cn    = nan(1, nt);
nfew  = 0;  nill = 0;

for it = 1:nt
    e  = eta(:,it);
    ok = isfinite(e);
    ngood(it) = nnz(ok);
    if ngood(it) < max(opt.minobs, nm+2)                  % need headroom over nm
        nfew = nfew + 1;  continue
    end

    Wk = sqrt(w(ok));                       % sqrt because LS minimises squares
    Aw = bsxfun(@times, A(ok,:), Wk);
    cn(it) = cond(Aw);
    if cn(it) > opt.maxcond                 % these levels cannot separate these
        nill = nill + 1;  continue          % modes -- see the note above
    end
    ew = e(ok).*Wk;
    ahat(:,it) = Aw\ew;
end

fprintf('eta fit: %d of %d timesteps solved (%d too few levels, %d ill-conditioned)\n', ...
        nnz(all(isfinite(ahat),1)), nt, nfew, nill);
if nill > 0
    fprintf('         rejected cond ranged %.3g to %.3g against a limit of %g;\n', ...
            min(cn(cn > opt.maxcond)), max(cn(cn > opt.maxcond)), opt.maxcond);
    fprintf('         median cond of the accepted fits was %.1f\n', median(cn(cn <= opt.maxcond),'omitnan'));
end

% --- cumulative reconstructions and diagnostics ----------------------------
recon = nan(nz, nt, nm);
for m = 1:nm
    recon(:,:,m) = A(:,1:m)*ahat(1:m,:);
end
resid  = eta - recon(:,:,nm);
varexp = nan(nm, nt);
for it = 1:nt
    e = eta(:,it); ok = isfinite(e);
    if nnz(ok) < max(opt.minobs, nm+2), continue; end
    v0 = sum(w(ok).*e(ok).^2);
    for m = 1:nm
        r = e(ok) - recon(ok,it,m);
        varexp(m,it) = 1 - sum(w(ok).*r.^2)/v0;
    end
end

fit.recon = recon;  fit.resid = resid;  fit.varexp = varexp;
fit.ngood = ngood;  fit.w = w;          fit.A = A;    fit.cond = cn;
end