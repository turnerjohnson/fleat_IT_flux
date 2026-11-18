% turner_tinker_M_DOD_moorings.m



clear; close all
addpath(genpath('/Users/tej/Documents'))

%% ADCP

M_ADCP = load("/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/FLEAT_M_ADCP_8m_BindDp.mat");

%% M1
M1_TSPbar = load("/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/M1_TSPbar.mat");


M1_TZdepth = load("/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/M1_TZdepth.mat");
 
M1.T = M1_TZdepth.TZdepth';
M1.time = M1_TZdepth.MTime;
M1.depth = M1_TZdepth.Zdepth;

%%
test_dod
% inside is dod_probe

%%
plot_dod_TS

%%

plot_adcp_time_misalign

%%

B = load("/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/M1_TSPbar.mat");
lat = 8.2;

M1.time  = B.MTime;                 % [1 x 44569]
M1.T     = B.TSdbar';               % [5 x 44569]
M1.S     = B.Sdbar';
M1.P     = B.PSdbar';               % dbar, raw
M1.depth = sw_dpth(M1.P, lat);      % m

figure(1);clf;
subplot(2,1,1)
psecline_mooring(M1, 'T', [], 'depth', [], struct('mooring','M1','lat',lat));

subplot(2,1,2)
psecline_mooring(M1, 'S', [], 'depth', [], struct('mooring','M1','lat',lat));

%%
M1_TZdepth = load("/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/M1_TZdepth.mat");
M1.T = M1_TZdepth.TZdepth';
M1.time = M1_TZdepth.MTime;
M1.depth = M1_TZdepth.Zdepth;
figure(2)
psectime_mooring(M1, 'T', '', [], [], [], 'depth', 1000)

%%

d = '/Users/tej/Documents/research_fleat/fleat_IT_flux/data/dod/';
m = 5;
A = load([d 'FLEAT_M_ADCP_8m_BindDp.mat'], sprintf('ADCPtm%d',m), ...
         sprintf('U%d',m), sprintf('V%d',m), 'ZZ');

M1a.time  = A.(sprintf('ADCPtm%d',m))(:).';          % [1 x 7547]
M1a.u     = A.(sprintf('U%d',m)) / 100;              % [75 x 7547] cm/s -> m/s
M1a.v     = A.(sprintf('V%d',m)) / 100;
M1a.depth = repmat(A.ZZ(:), 1, numel(M1a.time));     % ZZ = 8:8:600, fixed bins

cfg = struct('mooring', sprintf('M%d ADCP', m));
cl  = linspace(-0.6, 0.6, 15);                       % shared symmetric levels

figure(5); clf
subplot(2,1,1)
psectime_mooring(M1a, 'u', '', cl, [], [], 'depth', 1000, cfg);
subplot(2,1,2)
psectime_mooring(M1a, 'v', '', cl, [], [], 'depth', 1000, cfg);

%%
plot_adcp_aux