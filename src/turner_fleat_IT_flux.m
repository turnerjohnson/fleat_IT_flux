% turner_fleat_IT_flux.m
%
% THE ARRAY DRIVER. Internal tide energy flux <u'p'> = F at all five Palau
% FLEAT moorings, and the array-frame movie at the end.
%
% The methodology is in three layers and this file is only the top one:
%
%   fN_config         what is true of mooring N. Site, record window, T/S
%                     gates, mixed layer. Thin, because everything shared
%                     lives in fleat_config_base.
%   fleat_IT_flux     the pipeline, once, mooring-agnostic. cfg in, flux out.
%                     There is no mooring name anywhere inside it.
%   this file         call the configs, run the pipeline over the array,
%                     compare, animate.
%
% author TEJ, based on methods from Alford03 and Nash05, with help from Claude
% 15 October 2026

% --------------------------------------------
clear; close all
addpath(genpath('/Users/tej/Documents'))

%%
pmooring_layout
%h = pmooring_layout(struct('mode','gridded'));

%% ==== 0. what to run =================================================

names = {'F1','F2','F3','F4','F6'};     % the whole array ,'F2','F3', % 

cachedir  = '/Users/tej/Documents/research_fleat/fleat_IT_flux/data';
use_cache = false;     % reuse flux_<name>.mat if it is there

opt = struct();
opt.verbose   = true;
opt.keep      = 'flux';   % 'full' only when looking at one mooring, see cell 3
opt.showPlots = false;    % the per-mooring figures. Off for an array run.
opt.glance    = true;

%% ==== 1. run the array ===============================================
R = cell(1, numel(names));
for i = 1:numel(names)
    nm = names{i};
    cachefile = fullfile(cachedir, sprintf('flux_%s.mat', nm));

    if use_cache && isfile(cachefile)
        R{i} = load(cachefile).out;
        fprintf('[%s] from cache: %s\n', nm, cachefile);
        continue
    end

    % Each mooring's figures get their own hundred, so an array run with
    % showPlots on does not draw F6 over F1.
    o = opt;  o.figbase = 100*i;

    % catch for failed flux calc on one mooring.
    try
        cfg  = feval([lower(nm) '_config']);
        % big step here:
        R{i} = fleat_IT_flux(cfg, o);
        if ~isfolder(cachedir), mkdir(cachedir); end
        out = R{i};  save(cachefile, 'out');  clear out
    catch ME
        fprintf(2, '[%s] FAILED: %s\n', nm, ME.message);
        R{i} = struct('name', nm, 'ok', false, 'err', ME.message, ...
                      'lat', NaN, 'lon', NaN, 'H', NaN, 'time', [], 'flux', []);
    end
end

%% ==== 2. the array, side by side =====================================
% The comparison table. This is the point of the whole exercise: one pipeline,
% one set of assumptions, five sites, and the differences are the result.
%
% steadiness = |mean F| / mean |F|. Near 1 is a propagating wave with a fixed
% direction. Near 0 means the flux reverses and the record mean is hiding a
% standing wave. Alford finds O(1 kW/m) at his moorings -- compare.
fprintf('\n%-4s %8s %9s %6s %7s %7s %8s %6s %7s\n', ...
        'moor','lat','lon','H','|F|','dir','steady','c1','days');
fprintf('%s\n', repmat('-',1,70));
for i = 1:numel(R)
    r = R{i};
    if ~isstruct(r) || ~isfield(r,'ok') || ~r.ok
        fprintf('%-4s   (no flux)\n', r.name);  continue
    end
    d = r.diag;
    fprintf('%-4s %8.4f %9.4f %6.0f %7.3f %7.0f %8.2f %6.2f %7.0f\n', ...
            r.name, r.lat, r.lon, r.H, d.mag_mean/1e3, d.dir_mean, ...
            d.steadiness, d.c(1), d.ndays);
end
fprintf('\n|F| is mean magnitude in kW/m, dir is deg true, c1 in m/s.\n');

% Flag anything that ran on an unverified config, so the table above cannot be
% quoted without the caveat travelling with it.
for i = 1:numel(R)
    r = R{i};
    if isstruct(r) && isfield(r,'ok') && r.ok && ~isempty(r.cfg.unverified)
        fprintf('  %s ran with unverified: %s\n', r.name, strjoin(r.cfg.unverified, ', '));
    end
end

%% ==== 3. one mooring, with the diagnostics ===========================

dopt = opt;
dopt.keep        = 'full';    % attaches Mh, bkgd, eta, modes, ahat, uv
dopt.showPlots   = true;
dopt.figbase     = 0;
dopt.glance      = true;      % the T/S/u/v section overview
dopt.check_order = true;      % fit_first vs filter_first, by mode
% The sweep gate needs nmode_sweep+2 good levels in an hour, for every m, so it
% is capped by the LEVEL COUNT and not by curiosity. 10 against F6's 7 levels
% needed 12, no hour qualified, and every entry came back NaN. 4 needs 6, which
% F6 has. Raise it for a deeper mooring; F2's 40 levels will take 10 happily.
dopt.nmode_sweep = 4;         % true m-mode refits, raw and band
dopt.eta_spectra = true;      % peak above continuum, per mode
dopt.ellipse     = true;
dopt.mode_movie  = true;     % writes modes_<name>_M2.mp4
dopt.synth_date  = datenum(2016,10,15);   % synthetic NW arrival vs observed


cfg1 = f1_config();
cfg1.nbc = 3;              % fit_maxcond stays at the base 100
dopt.nmode_sweep = 8;
R3 = fleat_IT_flux(cfg1, dopt);

%%
puvp_slice(R2, datenum(2016,11,5,12,0,0));   % u', v', p', u'p'(z), and F
% puvp_slice(R2, 5000, 'obs', false);           % by index, no observed overlay
%
%
% switch N2 to gridded?:

% cfg = R2.cfg;  cfg.N2_source = 'gridded';
% R2b = fleat_IT_flux(cfg, dopt);
% fprintf('c1: sigma %.3f, gridded %.3f m/s\n', R2.diag.c(1), R2b.diag.c(1));

%% ==== 4. the array on one clock ======================================

A = align_flux_records(R, 1);

% Where the arrows point, on average, at each site. The static version of the
% movie, and the figure that should exist before the movie does.
figure('Color','w')
quiver(A.lon, A.lat, mean(A.Fx,2,'omitnan')', mean(A.Fy,2,'omitnan')', 0.6, ...
       'LineWidth', 1.6, 'MaxHeadSize', 0.6); hold on
plot(A.lon, A.lat, 'k.', 'MarkerSize', 14)
text(A.lon + 0.004, A.lat, A.name, 'FontSize', 9)
xlabel('longitude'); ylabel('latitude'); axis equal; box on; grid on
title('record-mean baroclinic M_2 energy flux')

%% ==== 5. <u'p'> AND <v'p'> THROUGH TIME ==============================
% The flux as a time series instead of a map arrow: two stacked panels, one
% per component -- <u'p'> on top, <v'p'> below -- with every mooring drawn in
% both and COLOUR CARRYING THE MOORING. A panel answers "do these sites
% agree?"; reading down the stack answers "which way is it going?". Both
% panels share one y-axis, since u and v are the same quantity on two axes
% and separate rulers would make a weak cross-shore flux look strong.
%
% The question it is for: is the baroclinic flux modulated at the spring/neap
% beat (14.765 d, M2 against S2)? Hence the 45-day default window -- three
% cycles, the shortest span in which a regular spacing is a claim and not an
% impression. That reasoning, and the warning that the dotted 14.765-d ruler's
% PHASE IS ARBITRARY (it is anchored at t0, not at syzygy -- nothing here knows
% the lunar phase), are in the header of plot/pflux_uvp_series.
%
% It takes R straight from section 1, or loads data/flux_<name>.mat itself in a
% fresh session. Both F2 and F3 run 2016-06-05 to 2017-04-13, so any start in
% there has full overlap -- but do not start ON the first day, where the
% bandpass and the 25-h boxcar both have edge effects. Sep-Nov 2016 sits well
% inside both records.

want = {'F2','F3', 'F4'};
t0   = datenum(2016,7,18);        % <-- the start to move around

if exist('R','var') && exist('names','var')
    isel = cellfun(@(w) find(strcmp(names, w), 1), want, 'UniformOutput', false);
    isel = cell2mat(isel(~cellfun(@isempty, isel)));
else
    isel = [];                   % nothing in memory -- load the cache
end

if isempty(isel)
    h5 = pflux_uvp_series(want, t0);
else
    h5 = pflux_uvp_series(R(isel), t0);
end

% ---- the knobs, all optional ----------------------------------------
% h5 = pflux_uvp_series(R(isel), t0, struct('ndays',90));            % 6 cycles
% h5 = pflux_uvp_series(R(isel), t0, struct('smooth_days',0, ...     % raw only
%                                           'ruler',false));
% h5 = pflux_uvp_series(R(isel), t0, struct('t_end',datenum(2016,12,1)));
% h5 = pflux_uvp_series(R(isel), t0, struct('sameY',false));         % per-panel y
% h5 = pflux_uvp_series(R(isel), t0, struct('raw',false));           % no dots,
%                                     % which also tightens the y-axis back onto
%                                     % the phase-averaged lines
% h5 = pflux_uvp_series(R, t0);                                      % all five

%% ==== 6. THE BAROTROPIC M2 FORCING, FROM TPXO ========================
% The same window, the same layout, the barotropic tide instead of the flux:
% pbt_tpxo_series draws u_BT on top and v_BT below out of
% data/u_m2_tpxo10_atlas_30_v2.nc, so figures 21 and 22 can be read against
% each other on screen.
%
% It needs nothing from the pipeline -- the configs give lat, lon and H -- so
% it runs in a cold session. TWO THINGS IT CANNOT SHOW, both in the header of
% core/tpxo_m2_bt and printed on the figure:
%
%   1. The atlas holds TRANSPORT (cm^2/s), so velocity is transport/H, and H
%      is the RECORDED depth because grid_tpxo10_atlas_30_v2.nc is not in
%      data/. That is a pure scale error on the amplitude: phase and shape are
%      untouched, so the series looks healthy while being tens of percent off.
%      Download the grid file and pass its h to fix it.
%
%   2. It is M2 ALONE, so there is NO SPRING/NEAP in it. The fortnight is the
%      M2-S2 beat, and S2 is a different file. The dashed envelope is flat on
%      purpose. So this figure can be compared with section 5 on PHASE, and
%      not on envelope, until u_s2_tpxo10_atlas_30_v2.nc is fetched and the
%      two predictions are summed.
%
% 45 days of a 12.42-h signal is 87 cycles and renders as a band. Drop to
% opt.ndays = 3 to see the waveform; keep 45 to sit beside the flux figure.

h6 = pbt_tpxo_series(want, t0);            % same want and t0 as section 5

% ---- the knobs ------------------------------------------------------
% h6 = pbt_tpxo_series(want, t0, struct('ndays',3));        % the waveform
% h6 = pbt_tpxo_series(want, t0, struct('H',[1515 1666]));  % explicit depths
% h6 = pbt_tpxo_series(want, t0, struct('file','/path/to/other_atlas.nc'));
%
% The node actually used, its distance from the mooring and the harmonic
% constants are on h6.B{i} and in the console summary. A node more than a few
% km away is flagged: at 1/30 deg the nearest WET node to a mooring on the
% Velasco flank need not be the mooring's own water column.

%% ==== 7. THE ARRAY MOVIE =============================================
% array_map_movie(R, mopt). Written 30 Oct 2026; the spec that used to live
% here is now the header of core/array_map_movie.m, decision by decision.
%
% It takes R and not A, because A carries only the flux. The velocity and the
% scalar dots come off out.uv, out.modes and out.Mh, which exist only under
% opt.keep = 'full' -- hence the two paths below.
%
% The scale, the missing-mooring rule, the 6-hour stride, the time-series
% cursor and the unverified-config footnote are all implemented as specified.
% Two things the spec did not anticipate, both in the function header:
%   - the bathymetry is banded on NON-UNIFORM levels and coloured by band
%     INDEX, because this box runs from reef flat to 3390 m and a uniform
%     contour renders the Velasco flank as one cliff;
%   - 'uv_mean' is the OBSERVED velocity, not u'. The depth mean of the
%     baroclinic modes is zero by construction, so a map of it is a map of
%     round-off.

% ---- the static frame first. This is the figure that has to exist. ----
S0 = array_map_movie(R, struct('box', [134.3 134.75 8.4 8.8]));

% ---- the movie: 6-hourly over the record, ~1250 frames ----
% mopt.outfile = '' animates on screen without writing the mp4.
mopt = struct();
mopt.box     = [134.3 134.75 8.4 8.8];
mopt.vector  = 'flux';
mopt.frames  = 'all';          % 'all' + stride_hr, or explicit indices
mopt.stride_hr = 6;
mopt.outfile = 'flux_array_M2.mp4';
mopt.fps     = 10;
mopt.scale   = 10 * 0.00144;

%%
M1 = array_map_movie(R, mopt);

%% ---- the other field choices ----------------------------------------
% All five need opt.keep = 'full' except 'flux' and 'flux_mean', so rerun
% section 1 with o.keep = 'full' before uncommenting these. They are much
% bigger in memory: Mh alone is the hourly field at every sensor.
%
% mopt.vector = 'flux_mean';                    % W/m^2, NOT the same ranking
% mopt.vector = 'bt';                           % barotropic M2 velocity
% mopt.vector = 'uv_slice';  mopt.z = 200;      % u'(200 m)
% mopt.vector = 'uv_mean';                      % observed, sampled-range mean
%
% ---- a scalar dot beside each mooring --------------------------------
% mopt.scalar   = 'T';     mopt.scalar_z = 200;   % nearest sensor to 200 m
% mopt.scalar   = 'eps';   mopt.scalar_z = [];    % sensor mean, auto log10
% mopt.scalar   = 'chi';
%
% chi and eps are log10-scaled automatically and T and S are not; override
% with mopt.scalar_log. The dots get their own colourbar under the map and
% are placed once, crowded moorings first, so they miss each other and the
% arrows -- F4 and F6 are 1.7 km apart and would otherwise overlap.
