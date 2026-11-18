function F_h_out = compute_density_mooring(F_h, cfg)
% COMPUTE_DENSITY_MOORING  Pressure and potential density on the mooring grid.
%
%   F_h = compute_density_mooring(F_h, cfg)
%
%   Adds .p (sw_pres from depth and cfg.lat) and .sigma (sw_pden referenced to
%   the surface, minus 1000).
%
%   ---------------------------------------------------------------------
%   WHY THERE IS A SALINITY CHECK BEFORE THE CALL
%   ---------------------------------------------------------------------
%   A NEGATIVE SALINITY MAKES sigma COMPLEX, and it does it silently. The chain
%   is sw_pden -> sw_dens -> sw_dens0 and sw_seck, and both of those take a
%   fractional power of S:
%
%       sw_dens0.m line 78    ... + (c0 + (c1 + c2*T68).*T68).*S.*sqrt(S) + ...
%       sw_seck.m  line 120   SR = sqrt(S);
%
%   sqrt of a negative number is imaginary, so one bad sample makes one element
%   of sigma complex -- and because a MATLAB array has a single type, THE WHOLE
%   ARRAY becomes complex and every good element gains a +0.000000i. That is why
%   the symptom looks global ("complex everywhere") when the cause is a handful
%   of samples: counting with nnz(imag(sigma) ~= 0) usually finds a few dozen
%   out of millions.
%
%   Nothing downstream complains either. sort, unique, min, max and the
%   relational operators all accept complex input, so the first real failure is
%   an interp1 or a plot several functions later, with no hint of where it came
%   from. F4's gridded product carries salinities from -114 to +489 psu at its
%   top five sensors; F2's never leave [34.33 35.07].
%
%   THE FIX IS qc_ts, NOT THIS FUNCTION. cfg.S_valid gates the raw samples
%   before the hourly average, which is the right place: one bad 10-minute value
%   inside an hourly mean contaminates the whole hour, and masked first the hour
%   is simply the mean of its five good neighbours. fleat_IT_flux already calls
%   qc_ts for every mooring. A script that calls this function directly without
%   having run qc_ts is the case this check exists to catch.
%
%   NOTE ON WHY THE HOURLY MEAN CAN HIDE IT: six 10-minute samples, five near
%   34.6 and one at -114, average to +9.8 -- positive, so sqrt is real and
%   sigma comes back real but WRONG, by about 19 kg/m^3. A complex sigma is the
%   loud version of this failure. The quiet version is worse, and only the gate
%   catches it.
%
%   The depth-side equivalent is ruled out rather than checked: sw_pres line 71
%   takes sqrt((1-C1).^2 - 8.84e-6*depth), which at this latitude only goes
%   negative for a depth beyond about 112 km.
%
% as part of pipeline for computing IT flux on FLEAT moorings
%
% author TEJ, with help from claude
% 29 September 2026, salinity guard 16 October 2026

%% ---- salinity has to be physical before it reaches sqrt -------------
S   = F_h.S;
neg = isfinite(S) & S < 0;
if any(neg(:))
    per  = sum(neg, 2);
    hit  = find(per > 0);
    zbar = mean(F_h.depth, 2, 'omitnan');
    msg  = sprintf(['%s: %d salinity sample(s) are NEGATIVE, which makes ' ...
                    'sw_pden return a COMPLEX sigma (sqrt(S) in sw_dens0 and ' ...
                    'sw_seck). Worst value %.2f psu.\n'], ...
                    mname(cfg), nnz(neg), min(S(neg)));
    for i = hit(:)'
        msg = [msg sprintf('    sensor %d (%.0f m): %d sample(s), min %.2f psu\n', ...
                           i, zbar(i), per(i), min(S(i,neg(i,:))))];
    end
    msg = [msg sprintf(['  Run qc_ts on the RAW samples before the hourly ' ...
                        'average:\n' ...
                        '      [F.T, F.S] = qc_ts(F.T, F.S, F.depth, cfg);\n' ...
                        '  cfg.S_valid is %s. fleat_IT_flux does this for ' ...
                        'you; a hand-rolled script has to.'], ...
                        rngstr(cfg, 'S_valid'))];
    error('compute_density_mooring:negS', '%s', msg);
end

%% ---- pressure and potential density ---------------------------------
F_h.p     = sw_pres(F_h.depth, cfg.lat);
F_h.sigma = sw_pden(F_h.S, F_h.T, F_h.p, 0) - 1000;

% Belt and braces. If sigma came back complex despite S being non-negative,
% something other than the salinity did it and the message should say so rather
% than let it travel. Nothing known produces this; it is here because the
% failure is invisible if it ever happens again for a new reason.
if ~isreal(F_h.sigma)
    nc = nnz(imag(F_h.sigma) ~= 0);
    error('compute_density_mooring:complex', ...
        ['sigma is complex at %d of %d elements but no salinity is negative. ' ...
         'Check cfg.lat (%.4f) and the depth range (%.1f to %.1f m) -- those ' ...
         'are the only other inputs.'], ...
         nc, numel(F_h.sigma), cfg.lat, ...
         min(F_h.depth(:)), max(F_h.depth(:)));
end

%% ---- how much of the density signal is salinity doing? --------------
sigma_fixedS = sw_pden(repmat(mean(F_h.S,2,'omitnan'),1,size(F_h.S,2)), ...
                       F_h.T, F_h.p, 0) - 1000;
fprintf('rms sigma difference, glider S vs time-mean S: %.4f kg m^-3\n', ...
        std(F_h.sigma(:) - sigma_fixedS(:), 'omitnan'));

F_h_out = F_h;
end

%% ======================================================================
function s = mname(cfg)
if isstruct(cfg) && isfield(cfg,'mooring') && ~isempty(cfg.mooring)
    s = cfg.mooring;
else
    s = 'mooring';
end
end

%% ======================================================================
function s = rngstr(cfg, f)
% What the config says the gate is, so the message names the actual value.
if ~isstruct(cfg) || ~isfield(cfg, f) || isempty(cfg.(f))
    s = sprintf('EMPTY (%s is not set, so qc_ts would do nothing)', f);
else
    s = mat2str(cfg.(f));
end
end
