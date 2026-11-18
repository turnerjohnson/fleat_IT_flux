function B = set_N2_source(B, src)
%SET_N2_SOURCE  Point B.N2full at one of the two stored N^2 profiles.
%
%   B = set_N2_source(B, 'sigma')     % record-mean S and T through sw_bfrq
%   B = set_N2_source(B, 'gridded')   % time median of the supplied F.N2 field
%
%   F_BACKGROUND_V2 always computes both and leaves them on the struct as
%   B.N2full_sigma and B.N2full_grid. This only decides which one sits in
%   B.N2full and therefore which one BUILD_MODES sees. Nothing is recomputed
%   and nothing else on B changes, which is what makes it safe to call
%   interactively after looking at the plot.
%
%   STALE AFTER THIS: anything already built from the previous N2full -- modes,
%   eigenspeeds, the eta fits, p', the fluxes. Rerun from the modes cell.

if ~isfield(B, 'N2full_sigma')
    error('set_N2_source:old', ...
        ['This background struct has no N2full_sigma, so it did not come from ' ...
         'f_background_v2. Rebuild it before switching sources.']);
end

switch lower(strtrim(src))
  case 'sigma'
    B.N2full = B.N2full_sigma;

  case 'gridded'
    if ~isfield(B,'N2full_grid') || isempty(B.N2full_grid)
        error('set_N2_source:nogrid', ...
            ['No gridded N^2 on this struct. f_background_v2 only builds it ' ...
             'when the mooring struct carries an N2 field of the same size as ' ...
             'depth. Check that F.N2 survived mooring_hourly_mean.']);
    end
    B.N2full = B.N2full_grid;

  otherwise
    error('set_N2_source:src', ...
        'N2 source must be ''sigma'' or ''gridded'', got ''%s''.', src);
end

B.N2_source = lower(strtrim(src));

% build_modes errors on this anyway, but it is cheaper to find out here.
if any(~isfinite(B.N2full)) || any(B.N2full <= 0)
    error('set_N2_source:bad', ...
        'The "%s" profile is not finite and positive everywhere.', B.N2_source);
end
end
