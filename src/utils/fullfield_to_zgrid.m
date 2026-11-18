function [Xg, info] = fullfield_to_zgrid(X, D, z, opt)
% FULLFIELD_TO_ZGRID  A knockdown-following field onto a fixed depth grid.
%
%   [Xg, info] = fullfield_to_zgrid(X, D, z)
%   [Xg, info] = fullfield_to_zgrid(X, D, z, opt)
%
%   X    [nl x nt]   the field, one row per LOGGER (not per depth)
%   D    [nl x nt]   that logger's INSTANTANEOUS depth, m positive down. This
%                    is fullfield.depth, and it moves: at F2 range(depth,2) is
%                    4.1 to 62.9 m of knockdown.
%   z    [1 x nz]    the fixed output grid, m positive down, ascending
%   opt  .maxgap     m. Do not interpolate across a vertical gap in the
%                    loggers wider than this. Default 1.5 * the median logger
%                    spacing in that column, COMPUTED per column, not typed.
%        .quiet      true to suppress the progress lines. Default false.
%
%   Xg   [nz x nt]   X on z, column by column, NO extrapolation past the
%                    shallowest or deepest logger that reported at that hour
%   info .cover      [nz x 1]  finite fraction per output level
%        .nskip      columns with fewer than 2 finite loggers
%        .maxgap     [1 x nt]  the gap limit used in each column
%        .ztop .zbot [1 x nt]  shallowest / deepest reporting logger
%
%   WHY THIS IS A LOOP OVER COLUMNS AND CANNOT BE ONE interp1 CALL
%   The source depths change every sample, so there is no shared breakpoint
%   vector. That is the whole difference between this product and the DOD one:
%   M*_TZdepth.mat has Zdepth = 60:10:1000 fixed, so build_dod_mooring gets
%   away with a single interp1. Here the grid IS the knockdown record, which
%   is exactly why the existing pipeline carries the knockdown NaN mask all
%   the way to the reconstructions -- and why moving to a fixed grid is the
%   restructuring this function exists to support.
%
%   WHY GAPS ARE NOT FILLED  (opt.maxgap)
%   The logger set is not uniformly spaced: at F2 it runs ~10 m apart through
%   the thermocline and then jumps (741 -> 889, 1289 -> 1389, 1429 -> ~1496).
%   Interpolating linearly onto a 4-10 m output grid across a 150 m hole
%   invents a dozen levels that look like measurements, carry equal weight in
%   a least-squares mode fit, and are a straight line drawn between two
%   unrelated depths. So the default refuses the hole: an output level whose
%   bracketing loggers are further apart than
%   1.5 * median(diff(depths that reported)) is left NaN. A level sitting
%   exactly ON a reporting logger survives regardless -- that one is a
%   measurement, not a fill. Set opt.maxgap = Inf to get the unguarded
%   behaviour and see what it was filling in.
%
%   THE TWO ENDS ARE NOT EXTRAPOLATED, EVER
%   Above the shallowest logger and below the deepest, Xg is NaN. This
%   function will not reach to the surface -- that is the modal fit's job, and
%   the point of the new structure is that the FIT, not the data, carries the
%   cap. Filling it here would hide the assumption inside what looks like data.
%
%   author TEJ, with help from Claude -- 17 November 2026

if nargin < 4 || ~isstruct(opt), opt = struct(); end
if ~isfield(opt,'maxgap') || isempty(opt.maxgap), opt.maxgap = []; end
if ~isfield(opt,'quiet')  || isempty(opt.quiet),  opt.quiet  = false; end

X = double(X);  D = double(D);
if ~isequal(size(X), size(D))
    error('fullfield_to_zgrid:shape', ...
        ['X is %s and D is %s. They must be the same size: every value needs ' ...
         'the depth it was measured at.'], mat2str(size(X)), mat2str(size(D)));
end
z  = double(z(:));
nz = numel(z);
nt = size(X, 2);
if nz > 1 && any(diff(z) <= 0)
    error('fullfield_to_zgrid:zorder', ...
        'z must be strictly ascending (m positive down); it is not.');
end
if any(D(:) < 0)
    error('fullfield_to_zgrid:dsign', ...
        ['D has %d negative values, so it is probably negative-down. This ' ...
         'function wants depth POSITIVE down, like z. Flip it before calling.'], ...
         nnz(D(:) < 0));
end

Xg          = nan(nz, nt);
info.nskip  = 0;
info.maxgap = nan(1, nt);
info.ztop   = nan(1, nt);
info.zbot   = nan(1, nt);

for k = 1:nt
    d = D(:,k);  x = X(:,k);
    g = isfinite(d) & isfinite(x);
    if sum(g) < 2, info.nskip = info.nskip + 1; continue; end

    [dk, is] = sort(d(g));
    xk = x(g);  xk = xk(is);
    [dk, iu] = unique(dk, 'first');         % two loggers at one depth: take one
    xk = xk(iu);
    if numel(dk) < 2, info.nskip = info.nskip + 1; continue; end

    info.ztop(k) = dk(1);
    info.zbot(k) = dk(end);

    % the gap limit for THIS column, from the loggers that actually reported
    if isempty(opt.maxgap)
        mg = 1.5 * median(diff(dk));
    else
        mg = opt.maxgap;
    end
    info.maxgap(k) = mg;

    col = interp1(dk, xk, z, 'linear', NaN);

    % Blank anything filled across a gap wider than mg. discretize gives the
    % bracketing interval of each output level (NaN outside the loggers), so
    % gaps(ib) is the distance the linear fill actually spanned.
    if isfinite(mg)
        gaps = diff(dk);
        ib   = discretize(z, dk);
        in   = isfinite(ib);
        wide = false(nz, 1);
        wide(in) = gaps(ib(in)) > mg;
        wide = wide & ~ismembertol(z, dk, 1e-9, 'DataScale', 1);
        col(wide) = NaN;
    end

    Xg(:,k) = col;
end

info.cover = mean(isfinite(Xg), 2);
if ~opt.quiet
    fprintf(['  fullfield_to_zgrid: %d loggers -> %d levels x %d columns, ' ...
             '%.1f%% finite, %d column(s) with <2 loggers\n'], ...
             size(X,1), nz, nt, 100*mean(isfinite(Xg(:))), info.nskip);
    fprintf(['    shallowest reporting logger %.1f-%.1f m, deepest ' ...
             '%.1f-%.1f m, gap limit %.1f-%.1f m\n'], ...
             min(info.ztop), max(info.ztop), min(info.zbot), max(info.zbot), ...
             min(info.maxgap), max(info.maxgap));
end
end
