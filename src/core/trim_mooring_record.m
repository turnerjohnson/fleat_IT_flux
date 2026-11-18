function F_out = trim_mooring_record(F, cfg)

% as part of pipeline for computing IT flux on FLEAT moorings.

% author TEJ, with help from claude
% 30 September 2026

%   F     RAW gridded mooring struct, straight off the load. Not F_h: this has
%         to run BEFORE mooring_hourly_mean, so that the background mean, N2,
%         the modes, eta and every fit are all built from the same window

%   cfg   from fX_config -- uses .t_start and .t_end, either of
%         which may be [] to leave that end alone.
%

t_start = [];  t_end = [];
if isfield(cfg,'t_start'), t_start = cfg.t_start; end
if isfield(cfg,'t_end'),   t_end   = cfg.t_end;   end

if ~isfield(F,'time') || isempty(F.time)
    error('trim_mooring_record:notime', ...
        'F.time is missing, so there is no clock to trim against.');
end

nt0 = numel(F.time);
t   = F.time(:).';
span = @(tt) sprintf('%s to %s (%.0f days, %d samples)', ...
        datestr(tt(1),'dd-mmm-yyyy'), datestr(tt(end),'dd-mmm-yyyy'), ...
        tt(end)-tt(1), numel(tt));

if isempty(t_start) && isempty(t_end)
    fprintf('record: %s, nothing trimmed (cfg.t_start and cfg.t_end are both [])\n', ...
            span(t));
    F_out = F;
    return
end

keept = true(1, nt0);
if ~isempty(t_start), keept = keept & t >= t_start; end
if ~isempty(t_end),   keept = keept & t <= t_end;   end

if ~any(keept)
    error('trim_mooring_record:empty', ...
        ['The window %s to %s keeps no samples. The record runs %s -- check ' ...
         'cfg.t_start and cfg.t_end are datenums on that scale.'], ...
        tstr(t_start,'start of record'), tstr(t_end,'end of record'), span(t));
end

if all(keept)
    fprintf('record: %s, nothing trimmed (the window contains it all)\n', span(t));
    F_out = F;
    return
end

% Second clock. ns0 must differ from nt0 or the column-count test below cannot
% tell them apart; at F2 and F4 it is 7512 against 45072.
keeps = [];  ns0 = NaN;
if isfield(F,'slopetime') && isnumeric(F.slopetime) && ~isempty(F.slopetime)
    ns0 = numel(F.slopetime);
    if ns0 == nt0
        warning('trim_mooring_record:sameclock', ...
            ['slopetime and time both have %d columns, so the two clocks ' ...
             'cannot be told apart by column count. slopetime-clock fields ' ...
             'are being cut with the sample-clock mask.'], nt0);
    else
        ts    = F.slopetime(:).';
        keeps = true(1, ns0);
        if ~isempty(t_start), keeps = keeps & ts >= t_start; end
        if ~isempty(t_end),   keeps = keeps & ts <= t_end;   end

        % The sample clock kept something, so the slope clock should too. If it
        % did not, the two do not overlap and cutting to nothing would be worse
        % than leaving them whole and saying so.
        if ~any(keeps)
            warning('trim_mooring_record:noslope', ...
                ['The window keeps %d of %d samples but none of the %d slope ' ...
                 'hours, so the two clocks do not overlap. slope* and slopetime ' ...
                 'are left untrimmed and will NOT line up with time.'], ...
                 nnz(keept), nt0, ns0);
            keeps = [];
        end
    end
end

% cut every field on either clock 
fn = fieldnames(F);
ncut_t = 0;  ncut_s = 0;
for i = 1:numel(fn)
    v = F.(fn{i});
    if size(v,2) == nt0
        F.(fn{i}) = v(:, keept);              % sample clock
        ncut_t = ncut_t + 1;
    elseif ~isempty(keeps) && size(v,2) == ns0
        F.(fn{i}) = v(:, keeps);              % slopetime clock
        ncut_s = ncut_s + 1;
    end                                       % else static, left alone
end

% report 
fprintf('record: %s, trimmed from %d samples (%d fields on the sample clock', ...
        span(F.time(:).'), nt0, ncut_t);
if ncut_s > 0
    fprintf(', %d on the slope clock -> %d hours', ncut_s, nnz(keeps));
elseif isfinite(ns0)
    fprintf(', slope clock NOT cut');
end
fprintf(')\n');

F_out = F;
end

%% ======================================================================
function s = tstr(x, fallback)
% A datenum for the message, or the words for the end that was left open.
if isempty(x), s = fallback; else, s = datestr(x,'dd-mmm-yyyy'); end
end
