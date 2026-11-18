function [sig, info] = stabilizeSigma(sig, w, minstep)
%STABILIZESIGMA  Nearest stable (non-decreasing) density profile.
%   SIG = STABILIZESIGMA(SIG) treats each column of SIG (nz x n, ordered
%   shallow to deep) as a profile and returns the closest profile that is
%   monotonically increasing with depth, in a weighted least-squares sense.
%   Solved by pool-adjacent-violators; equivalent to the lsqlin formulation
%   of Rudnick's computeStableDensity but exact, toolbox-free and O(nz).
%
%   SIG = STABILIZESIGMA(SIG,W) weights level k by W(k) (e.g. bin counts),
%   so poorly-sampled levels take up more of the adjustment. Default: ones.
%
%   SIG = STABILIZESIGMA(SIG,W,MINSTEP) additionally forces a strict
%   increase of at least MINSTEP between levels. Use sparingly: flat
%   segments are a real statement that stratification is unresolved, and
%   MINSTEP invents structure. Default 0.
%
%   INFO returns .nfix (columns altered), .maxadj (largest change per
%   column) and .rmsadj, which are the numbers you should actually look at.
%
%   NaNs are handled per column: finite entries are stabilized among
%   themselves, NaNs pass through untouched.

if nargin < 2 || isempty(w),       w = ones(size(sig,1),1);  end
if nargin < 3 || isempty(minstep), minstep = 0;              end

[nz, n]   = size(sig);
w         = w(:);  w(~isfinite(w) | w <= 0) = eps;
info.nfix = 0;  info.maxadj = zeros(1,n);  info.rmsadj = zeros(1,n);

for j = 1:n
    y  = sig(:,j);
    ok = isfinite(y);
    if nnz(ok) < 2, continue; end
    yo = y(ok);  wo = w(ok);

    if minstep > 0                       % shift to an equivalent flat problem
        off = minstep * (0:nnz(ok)-1)';
        yo  = yo - off;
    end
    if all(diff(yo) >= 0), continue; end  % already stable

    x = pava(yo, wo);
    if minstep > 0, x = x + off; end

    d = x - y(ok);
    info.nfix = info.nfix + 1;
    info.maxadj(j) = max(abs(d));
    info.rmsadj(j) = sqrt(mean(d.^2));
    y(ok) = x;  sig(:,j) = y;
end
end

function x = pava(y, w)
% Pool adjacent violators: weighted non-decreasing least-squares fit.
n = numel(y);
val = zeros(n,1); wt = zeros(n,1); en = zeros(n,1);   % block value, weight, end index
j = 0;
for i = 1:n
    j = j + 1;  val(j) = y(i);  wt(j) = w(i);  en(j) = i;
    while j > 1 && val(j-1) > val(j)                   % merge with block above
        val(j-1) = (wt(j-1)*val(j-1) + wt(j)*val(j)) / (wt(j-1) + wt(j));
        wt(j-1)  =  wt(j-1) + wt(j);
        en(j-1)  =  en(j);
        j = j - 1;
    end
end
x = zeros(n,1);  lo = 1;
for b = 1:j
    x(lo:en(b)) = val(b);  lo = en(b) + 1;
end
end