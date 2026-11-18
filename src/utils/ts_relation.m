function [tc, sc, info] = ts_relation(t, S, opt)
% TS_RELATION  Binned S(T) curve for the DOD moorings, with its gaps marked.
%
%   [tc, sc]       = ts_relation(t, S)
%   [tc, sc, info] = ts_relation(t, S, opt)
%
%   t, S  pooled in-situ temperature and practical salinity from the
%         instruments. Any shape; NaNs ignored.
%   opt  .nbin     bins. Default 60.
%        .minn     minimum points to keep a bin. Default 50.
%        .smooth   moving-median then moving-mean window, in bins. Default 5.
%                  0 disables smoothing.
%        .nfine    fine equal-width bins used to FIND the voids. Default 400.
%        .gapcount a fine bin with this many points or fewer counts as empty.
%                  Default 0, i.e. strictly no water at that temperature.
%        .gapmin   shortest run of empty fine bins that counts as a void,
%                  in degC. Default 0.2.
%
%   tc, sc  the curve, ascending in temperature.
%   info .n          points per bin
%        .width      temperature width of each bin
%        .gap        logical, bins that span a void
%        .gap_range  [lo hi] temperature intervals with no real support
%
%   WHY EQUAL-COUNT BINS AND NOT EQUAL-WIDTH
%   This is the whole fix. The temperature distribution at these moorings is
%   not remotely uniform: the deepest sensor parks 44% of all samples inside
%   1.35 degC while the thermocline sensors spread the rest over 20 degC. With
%   80 equal-width bins the occupancy ran from 0 to 89,014 points, so the
%   median was rock-solid in some bins and built from a handful of samples in
%   others, and the curve inherited exactly that -- smooth where the data was
%   dense, visibly jagged between 11 and 13 degC where each bin held a few
%   thousand instead of tens of thousands.
%
%   Equal-count bins put the same number of samples behind every median, so
%   the curve has uniform reliability along its length. The cost is that bins
%   are no longer evenly spaced in temperature, which is the point: where the
%   data is sparse the bin simply gets wider.
%
%   THE VOID BETWEEN THE SENSORS IS NOT SMOOTHED OVER, IT IS MARKED
%   Pooled over all five
%   moorings there is no water at all sampled between 7.06 and 7.87 degC --
%   that is the 290-720 dbar hole between the 4th and deepest instruments,
%   500 dbar of it at M1. In the old equal-width binning that void showed up
%   as two under-occupied bins that were silently dropped, and interp1 then
%   drew a straight line across it. Here it is found by asking a fine
%   histogram which temperatures hold no water, and reported in .gap_range.
%
%   That matters past cosmetics. s_method 'ts' in build_dod_mooring maps this
%   curve through the gridded T field, so any water in the void would be
%   assigned a salinity read off a straight line drawn between two unrelated
%   water masses -- a fabricated salinity gradient at exactly the depths where
%   N2 feeds the mode fit. The curve is still returned across the gap, because
%   interpolating something is often the practical choice, but the caller is
%   told where it is unsupported and can NaN it instead. Nothing here decides
%   that silently.
%
%   SMOOTHING IS APPLIED AFTER BINNING, NOT INSTEAD OF IT
%   A moving median then a moving mean over the bin medians. The moving median
%   kills single-bin outliers without rounding off real curvature; the moving
%   mean then takes the staircase out. Both run over BINS, so with equal-count
%   bins they smooth over a constant number of samples rather than a constant
%   temperature interval.

%% ---- arguments ------------------------------------------------------
if nargin < 3 || ~isstruct(opt), opt = struct(); end
if ~isfield(opt,'nbin')   || isempty(opt.nbin),   opt.nbin   = 60; end
if ~isfield(opt,'minn')   || isempty(opt.minn),   opt.minn   = 50; end
if ~isfield(opt,'smooth') || isempty(opt.smooth), opt.smooth = 5;  end
if ~isfield(opt,'nfine')    || isempty(opt.nfine),    opt.nfine    = 400; end
if ~isfield(opt,'gapcount') || isempty(opt.gapcount), opt.gapcount = 0;   end
if ~isfield(opt,'gapmin')   || isempty(opt.gapmin),   opt.gapmin   = 0.2; end

t = double(t(:));  S = double(S(:));
g = isfinite(t) & isfinite(S);
t = t(g);  S = S(g);
if numel(t) < 10*opt.nbin
    error('ts_relation:few', ...
        'Only %d finite T/S pairs for %d bins.', numel(t), opt.nbin);
end

%% ---- equal-count bin edges -----------------------------------------
% By sorting rather than quantile(), so this does not need the Statistics
% toolbox. Duplicate edges are collapsed, which happens when one temperature
% is so over-represented that it spans several bins' worth of samples.
ts   = sort(t);
idx  = unique(round(linspace(1, numel(ts), opt.nbin+1)));
edges = unique(ts(idx)).';
if numel(edges) < 4
    error('ts_relation:degenerate', ...
        'Temperature takes too few distinct values to bin (%d edges).', numel(edges));
end
edges(1)   = edges(1)   - eps(edges(1))*10;        % keep the end points in
edges(end) = edges(end) + eps(edges(end))*10;

ib = discretize(t, edges);
nb = numel(edges) - 1;
[tc, sc, info.n, info.width] = deal(nan(1, nb));
for b = 1:nb
    k = ib == b;
    info.n(b)     = sum(k);
    info.width(b) = edges(b+1) - edges(b);
    if info.n(b) >= opt.minn
        tc(b) = median(t(k));
        sc(b) = median(S(k));
    end
end

%% ---- mark the voids -------------------------------------------------
% Asked directly of the data, on a FINE EQUAL-WIDTH histogram, rather than
% inferred from the equal-count bin widths. Those widths legitimately vary a
% hundredfold here -- 0.06 degC inside the deep blob to 2.24 degC across the
% thermocline's thin tail -- so "unusually wide bin" flags sparse-but-real
% stretches as readily as true voids. Asking which temperatures have NO water
% at all is the question actually worth answering, and it has one answer here:
% 7.06-7.87 degC, thirteen consecutive empty bins.
fe = linspace(min(t), max(t), opt.nfine + 1);
ce = histcounts(t, fe);
empty = ce <= opt.gapcount;
info.gap_range = zeros(0, 2);
b = 1;
while b <= numel(empty)
    if empty(b)
        e2 = b;
        while e2 < numel(empty) && empty(e2+1), e2 = e2 + 1; end
        if fe(e2+1) - fe(b) >= opt.gapmin
            info.gap_range(end+1, :) = [fe(b) fe(e2+1)];   %#ok<AGROW>
        end
        b = e2 + 1;
    else
        b = b + 1;
    end
end
info.gap = false(1, nb);
for r = 1:size(info.gap_range, 1)
    info.gap = info.gap | (edges(1:end-1) < info.gap_range(r,2) & ...
                           edges(2:end)   > info.gap_range(r,1));
end

%% ---- keep, smooth, return ------------------------------------------
keep = isfinite(tc) & isfinite(sc);
tc = tc(keep);  sc = sc(keep);
info.n = info.n(keep);  info.width = info.width(keep);  info.gap = info.gap(keep);
[tc, k] = sort(tc);  sc = sc(k);  info.gap = info.gap(k);
info.n = info.n(k);   info.width = info.width(k);

if opt.smooth > 1 && numel(sc) > opt.smooth
    sc = smoothdata(sc, 'movmedian', opt.smooth);
    sc = smoothdata(sc, 'movmean',   opt.smooth);
end
info.nbin_used = numel(tc);
end
