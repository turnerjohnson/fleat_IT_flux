function [F_h, info] = mooring_hourly_mean(F, cfg)
% MOORING_HOURLY_MEAN  Hourly block mean of every time-varying field.
%
%   F_h          = mooring_hourly_mean(F, cfg)
%   [F_h, info]  = mooring_hourly_mean(F, cfg)
%
%   F     gridded mooring data (Anna Simpson's fullfield from FLEAT), or the
%         synthetic equivalent from synth_baroclinic.
%   cfg   from f2_config or f4_config -- uses .avg_k (samples per hour).
%
%   F_h   the same struct on the hourly clock. Every field keeps its name and
%         its position in the struct, so F_h is a drop-in replacement for F in
%         psectime_mooring and psecline_mooring.
%   info  .averaged .passed_hourly .passed_static  cellstr of field names
%         .nh        hours out
%         .ndropped  samples lost to the partial final hour
%
%   WHY THIS IS GENERIC
%   The cell this replaces hand-listed six fields -- time, depth, T, S, u, v --
%   and silently dropped the other ten that are on the sample clock: dTdz, N2,
%   spd, chi, chiSSF, chi_nomask, eps, epsSSF, eps_nomask and slope20020.
%   Anything that wanted hourly N2 had to reach back to the raw struct and
%   re-average by hand, which is how the drivers ended up plotting a 10-minute
%   median against an hourly depth.
%
%   WHICH FIELDS GET AVERAGED
%   Column count against numel(F.time), and nothing else. That is the same test
%   pick_time uses in psectime_mooring and psecline_mooring, and the same one
%   the record trim in turner_fleat4_IT_flux uses.
%
%   Name-based rules do not work here. slope20020 is [nsens x 45072], on the
%   10-MINUTE clock, despite looking like one of the hourly slope fields; and
%   slopetime would itself match a 'slope*' prefix. The column count gets both
%   right without knowing anything about the product.
%
%   id and stats are [1 x nsens] cells -- transposed relative to every data
%   field -- so a rule keyed on "rows == nsens" would misread them as data. The
%   column count does not.
%
%   THE TWO CLOCKS LAND ON TOP OF EACH OTHER
%   45072/6 = 7512 = numel(slopetime), and slopetime agrees with the hourly mean
%   of time to better than a thousandth of a minute in both F2 and F4. So the
%   slope fields, which the product already computed hourly, arrive on exactly
%   the grid this function produces and are passed through untouched.
%
%   That stops being true after a record trim. cfg.t_end at F4 cuts time to 164
%   of 313 days, and the trim loop only subsets fields on the sample clock, so
%   slopetime is left spanning the whole deployment. Passing it through then
%   would hand back an F_h whose slope columns do not line up with its own
%   time. Instead the slopetime-clock fields are subset to the hours the output
%   actually spans. The values are unchanged; only the matching columns survive.
%
%   ONE MEAN FOR EVERYTHING
%   Arithmetic, over the k samples in each hour, NaNs omitted. Right for chi and
%   eps, which are fluxes. Conventional for N2. If you ever decide a log-normal
%   field wants a geometric mean in time, that is a deliberate per-field choice
%   to add here -- the report below names all sixteen averaged fields so the
%   assumption is visible rather than implicit.

%% ---- arguments ------------------------------------------------------
if nargin < 2 || ~isstruct(cfg)
    error('mooring_hourly_mean:nocfg', ...
        'Needs a config: mooring_hourly_mean(F, cfg), with cfg.avg_k set.');
end
if ~isfield(cfg,'avg_k') || ~isscalar(cfg.avg_k) || ~isfinite(cfg.avg_k) ...
        || cfg.avg_k < 1 || mod(cfg.avg_k,1) ~= 0
    error('mooring_hourly_mean:avgk', ...
        'cfg.avg_k must be a positive whole number of samples per hour.');
end
if ~isfield(F,'time') || isempty(F.time)
    error('mooring_hourly_mean:notime', ...
        ['F.time is missing, so there is no clock to test the other fields ' ...
         'against and nothing can be classified.']);
end

k  = cfg.avg_k;
nt = numel(F.time);
nh = floor(nt/k);
if nh < 1
    error('mooring_hourly_mean:short', ...
        'Only %d samples, which is less than the %d in one hour.', nt, k);
end
info.nh       = nh;
info.ndropped = nt - nh*k;

% Straight to [nrow x nh]. The version of this in the drivers squeezed first,
% which collapses a [1 x nt] field such as time to a column, and relied on a
% trailing reshape to put it back; this needs no undoing and is also correct
% when there is only one row.
blockmean = @(X) reshape(mean(reshape(X(:,1:nh*k), size(X,1), k, nh), ...
                              2, 'omitnan'), size(X,1), nh);

%% ---- classify and average every field -------------------------------
[info.averaged, info.passed_hourly, info.passed_static] = deal({});

fn = fieldnames(F);
for i = 1:numel(fn)
    f = fn{i};
    v = F.(f);

    if (isnumeric(v) || islogical(v)) && ndims(v) == 2 && size(v,2) == nt
        F_h.(f) = blockmean(double(v));         % double() so logicals average
        info.averaged{end+1} = f;               % to a fraction of the hour
    else
        F_h.(f) = v;
        if isnumeric(v) && ndims(v) == 2 && size(v,2) ~= 1 && size(v,2) ~= nt
            info.passed_hourly{end+1} = f;      % candidate second clock
        else
            info.passed_static{end+1} = f;
        end
    end
end

if ~ismember('time', info.averaged)
    error('mooring_hourly_mean:timeclass', ...
        ['F.time was not classified as a sample-clock field. It must be a ' ...
         'numeric [1 x nt] row; got %s %s.'], class(F.time), mat2str(size(F.time)));
end
t_h = F_h.time;

%% ---- keep the second clock aligned ----------------------------------
% Only fires when the record was trimmed. On a full F2 or F4 record the slope
% fields already have nh columns and drop straight through.
if ~isempty(info.passed_hourly)
    tol  = 0.5/24;                                  % half an hour, in days
    keep = [];
    if isfield(F,'slopetime') && isnumeric(F.slopetime)
        keep = F.slopetime(:).' >= t_h(1) - tol & F.slopetime(:).' <= t_h(end) + tol;
    end

    for i = 1:numel(info.passed_hourly)
        f  = info.passed_hourly{i};
        nc = size(F.(f), 2);
        if nc == nh, continue; end                  % already on the output grid

        if ~isempty(keep) && numel(keep) == nc && nnz(keep) == nh
            F_h.(f) = F.(f)(:, keep);
            fprintf(['mooring_hourly_mean: %s subset from %d to %d columns to ' ...
                     'match the trimmed record\n'], f, nc, nh);
        else
            warning('mooring_hourly_mean:clock', ...
                ['%s has %d columns, which is neither the %d samples of F.time ' ...
                 'nor the %d hours out. It is passed through unchanged and does ' ...
                 'NOT line up with F_h.time -- do not plot the two against each ' ...
                 'other without interpolating first.'], f, nc, nt, nh);
        end
    end
end

%% ---- report ---------------------------------------------------------
fprintf('hourly: %d sensors, %d hours (%.0f days), %d samples dropped\n', ...
        nsens_of(F_h), nh, t_h(end)-t_h(1), info.ndropped);
report('averaged',      info.averaged);
report('hourly, passed', info.passed_hourly);
report('static, passed', info.passed_static);
end

%% ======================================================================
function report(label, names)
if isempty(names), return; end
fprintf('  %s (%d): %s\n', label, numel(names), strjoin(names, ' '));
end

%% ======================================================================
function n = nsens_of(F_h)
% Rows of the first averaged field that has more than one, so the count means
% sensors rather than the single row of time.
n = 1;
fn = fieldnames(F_h);
for i = 1:numel(fn)
    v = F_h.(fn{i});
    if isnumeric(v) && size(v,1) > 1, n = size(v,1); return; end
end
end
