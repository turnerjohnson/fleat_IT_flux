function [ys, isfull] = runmean_z(y, z, L)
%RUNMEAN_Z  Running mean over a window L metres wide, on an irregular grid.
%
%   This is comparable to nsqSolo's
%
%       ii = floor(nzsmooth/2) + 1;  jj = mz - floor(nzsmooth/2);
%       Nsqs(ii:jj,:) = nanrunmean(bindata.Nsq, nzsmooth);
%
%   written for levels that are not evenly spaced. For a uniform grid of
%   spacing dz and L = (n-1)*dz with n odd the two are IDENTICAL: the same
%   points fall in each window, the mean is unweighted the same way, and
%   ISFULL is true on exactly ii:jj.
%
%   The difference is what happens when the spacing varies. A fixed point
%   count is a different physical filter at every depth -- on F2, 5 points is
%   50 m at 200 m and 500 m at 800 m. A fixed width is the same filter
%   everywhere, and in a sparse stretch its window simply holds one point and
%   returns it unchanged. That is the honest answer: you cannot smooth
%   structure you never sampled.
%
%   ISFULL marks levels whose window fits entirely inside the data. Those are the
%   ones nsqSolo keeps. The rest are returned UNSMOOTHED rather than NaN,
%   because build_modes needs a finite column, and a shrinking window would
%   quietly make the ends a different estimator from the middle.

% October 2026, TEJ w/ help from claude

z  = z(:);  y = y(:);
n  = numel(z);
ys     = y;
isfull = false(n,1);
h  = L/2;
tol = 1e-9;

for i = 1:n
    isfull(i) = (z(i) - h >= z(1) - tol) && (z(i) + h <= z(end) + tol);
    if ~isfull(i), continue; end
    w = abs(z - z(i)) <= h + tol;
    if any(isfinite(y(w)))
        ys(i) = mean(y(w), 'omitnan');
    end
end
end