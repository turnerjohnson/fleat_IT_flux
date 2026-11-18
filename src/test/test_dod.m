% test_dod.m -- smoke test for dod_probe + build_dod_mooring
clear; close all
addpath(genpath('/Users/tej/Documents/research_fleat/fleat_IT_flux/src'))
addpath('/Users/tej/Documents/MATLAB/gsw')

fprintf('\n########## 1. PROBE ##########\n');
info = dod_probe;

fprintf('\n########## 2. BUILD M1 ##########\n');
[F, bi] = build_dod_mooring(1);
whosstruct(F)   % no semicolon on purpose: this is the field/size table
fprintf('\nreadme:\n');  fprintf('  %s\n', F.readme{:});

fprintf('\n########## 3. ASSERTIONS ##########\n');
res = {'FAIL','PASS'};
% all(c(:)) so a condition that comes back non-scalar collapses to one
% logical instead of indexing res with a vector -- that returns a
% comma-separated list and fprintf then recycles its format over every
% element, which looks like runaway output.
chk = @(c, s) fprintf('  [%s] %s\n', res{double(all(logical(c(:))))+1}, s);
chk(abs(median(diff(F.time))*24 - 1) < 1e-9,            'time step is exactly 1 hour');
chk(all(abs(mod(F.time*24, 1)) < 1e-9),                 'time lands on exact hours');
chk(size(F.T,2) == numel(F.time),                       'T is [nz x nt]  (not transposed)');
chk(size(F.T,1) == numel(F.z),                          'T rows == numel(z)');
chk(isequal(size(F.u), size(F.T)),                      'u matches T in shape');
chk(max(abs(F.u(:))) < 2,                               'u looks like m/s not cm/s');
chk(mean(isfinite(F.T(:))) > 0.9,                       'T coverage > 90%');
chk(mean(isfinite(F.u(:))) > 0.8,                       'u coverage > 80%');
chk(median(F.N2(isfinite(F.N2))) > 1e-6 && ...
    median(F.N2(isfinite(F.N2))) < 1e-3,                'N2 magnitude is plausible');
chk(median(F.T(1,:),'omitnan') > median(F.T(end,:),'omitnan'), 'T decreases with depth');
chk(size(F.inst.S,1) == bi.nsens,                       'inst rows == instrument count read from file');
chk(ismember(bi.nsens, [5 6]),                          'instrument count is 5 or 6');
chk(max(F.S(isfinite(F.S))) < 35.6,                     'no S spikes survived qc');
chk(min(F.inst.C(isfinite(F.inst.C))) > 2.5,            'no negative conductivity survived qc');
chk(isequal(F.z, 10:10:1000),                           'grid is the union range 10:10:1000');
chk(all(isnan(F.T(1,:))),                               'T is NaN at 10 m (no sensor reaches there)');
chk(mean(isfinite(F.u(1,:))) > 0.9,                     'u IS present at 10 m (the point of a union grid)');
chk(all(isnan(F.u(F.z > 600, :)), 'all'),               'u is NaN below the 600 m ADCP limit');
chk(mean(isfinite(F.T(F.z == 1000, :))) > 0.95,         'T survives to 1000 m');
chk(mean(isfinite(F.N2(F.z == 1000, :))) > 0.9,         'N2 survives to 1000 m via the deep-S extension');
for f = {'T','S','u','v','N2'}
    fprintf('  usable range %-3s: %4.0f - %4.0f m\n', f{1}, ...
        bi.range.(f{1})(1), bi.range.(f{1})(2));
end
R = info.identity;
chk(all(diag(R)' == max(R, [], 2)'),                    'ADCP suffix m matches T/S mooring m at all five');
chk(min(diag(R)) > 0.95,                                'all five knockdown diagonals > 0.95');

fprintf('\n########## 4. THE TWO uv METHODS, COMPARED ##########\n');
Fi = build_dod_mooring(1, struct('uv_method','interp'));
Fb = build_dod_mooring(1, struct('uv_method','binmean'));
d  = Fi.u - Fb.u;
g  = isfinite(d);
fprintf('  rms(interp - binmean) = %.4f m/s, = %.1f%% of rms(u), over %d cells\n', ...
    std(d(g)), 100*std(d(g))/std(Fi.u(isfinite(Fi.u))), sum(g(:)));

fprintf('\n########## 5. THE TWO S METHODS, COMPARED ##########\n');
Fp = build_dod_mooring(1, struct('s_method','pressure'));
Ft = build_dod_mooring(1, struct('s_method','ts'));
assert(isequal(size(Fp.S), size(Ft.S)), 'the two S fields are not the same shape');
d  = Fp.S - Ft.S;
g  = isfinite(d);                       % finite in BOTH methods
fprintf('  rms(pressure - ts) = %.4f psu over %.1f%% common coverage (%d of %d cells)\n', ...
    std(d(g)), 100*mean(g(:)), sum(g(:)), numel(g));
fprintf('  pressure-only cells %d, ts-only cells %d\n', ...
    sum(isfinite(Fp.S(:)) & ~isfinite(Ft.S(:))), ...
    sum(isfinite(Ft.S(:)) & ~isfinite(Fp.S(:))));

fprintf('\n########## 6. ALL FIVE ##########\n');
fprintf('  (watch the instrument counts: M2 and M3 have 6, the others 5)\n');
for m = 1:5
    [G.(sprintf('M%d', m)), gi(m)] = build_dod_mooring(m);
end
fprintf('\n  instrument counts: %s\n', mat2str([gi.nsens]));
fprintf('\ndone. G.M1 .. G.M5 built.\n');

fprintf('\n########## 7. EYEBALL IT ##########\n');
figure('Position',[100 100 1100 800]);
subplot(4,1,1); pcolor(F.time, F.z, F.T);  shading flat; axis ij
  datetick('x','mmm'); ylabel('z (m)'); title('M1 T (degC)'); colorbar
subplot(4,1,2); pcolor(F.time, F.z, F.S);  shading flat; axis ij
  datetick('x','mmm'); ylabel('z (m)'); title('M1 S (psu)'); colorbar
subplot(4,1,3); pcolor(F.time, F.z, F.u);  shading flat; axis ij
  datetick('x','mmm'); ylabel('z (m)'); title('M1 u (m/s)'); colorbar
  caxis([-0.5 0.5])
subplot(4,1,4); plot(F.time, F.adcp.depth, 'k', F.time, F.adcp.Tdn, 'r')
  datetick('x','mmm'); legend('adcp.depth (m)','adcp.Tdn (degC, the "tilt" field)')
  title('knockdown and the mislabelled temperature')
