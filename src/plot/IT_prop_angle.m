% IT_prop_angle.m

%% --- M2 internal tide seen at one mooring depth ---
lat  = 7.5;                             % deg N (Palau)
f    = 2*(7.2921e-5)*sind(lat);         % Coriolis
om   = 2*pi/(12.4206*3600);             % M2 angular frequency
A    = 215;                               % propagation azimuth (toward), deg CW from NORTH
                                        %   0 = toward N, 90 = toward E, 135 = toward SE
                                        %   wave COMES FROM A+180
% A is TOWARD, the same convention as synth_baroclinic's P.to_deg,
% compute_tidal_flux's F.dir and compare_synth_real_uvp's to_deg. This wave is
% compare_synth_real_uvp's synthetic for a single mode with psi = pi/2: there
% th = -om*t + psi, so sin(th) = cos(om*t), the upar used here.
U0   = 0.05;                            % along-wave velocity amplitude [m/s]
rho0 = 1025;                            % reference density [kg/m^3]
c1   = 2.5;                             % mode-1 eigenspeed [m/s], sets the impedance
t    = linspace(0, 2*12.4206*3600, 500);% two M2 periods [s]

% clean polarization: along-wave in phase w/ p, cross-wave 90 deg & f/om smaller
upar = U0*cos(om*t);                    % along propagation
uprp = (f/om)*U0*sin(om*t);             % cross (rotary) part

% pressure. For a free mode the momentum equations give
%   upar = om*k*p'/(rho0*(om^2-f^2)),  with  om^2-f^2 = c1^2 k^2,
% so p' is IN PHASE with upar and the amplitudes are tied:
P0 = rho0*c1*U0*sqrt(1 - (f/om)^2);     % [Pa]
p  = P0*cos(om*t);

% rotate into geographic east (u) / north (v). Propagation points along
% (sind A, cosd A); the cross-wave axis is 90 deg CLOCKWISE of that,
% (cosd A, -sind A), which is what makes the vector turn clockwise -- the
% correct sense for the northern hemisphere.
u =  upar*sind(A) + uprp*cosd(A);
v =  upar*cosd(A) - uprp*sind(A);

% in phase => nonzero mean flux, directed along A
fprintf('p0 = %.1f Pa for u0 = %.3f m/s  (p0/u0 = rho0*c1*sqrt(1-(f/om)^2) = %.0f Pa per m/s)\n', ...
        P0, U0, P0/U0);
fprintf('<upar*p''> = %.3f W/m^2 toward %g deg CW from north\n', 0.5*U0*P0, A);

figure;
subplot(1,2,1)                          % propagation direction in lon/lat plane
quiver(0,0, sind(A), cosd(A), 0,'k','LineWidth',2,'MaxHeadSize',1); hold on
axis equal; axis([-1 1 -1 1]); grid on
xlabel('lon'); ylabel('lat');
title(sprintf('wave toward %g\\circ',A))

subplot(1,2,2)                          % what the mooring records
yyaxis left                             % velocity and pressure differ by O(10^3),
h1 = plot(t/3600,u,'-','LineWidth',1.5,'Color',[0.850 0.325 0.098]); hold on
h2 = plot(t/3600,v,'-','LineWidth',1.5,'Color',[0 0.447 0.741]);
ylabel('velocity (m s^{-1})')
yyaxis right                            %   so p' gets its own axis
h3 = plot(t/3600,p,'-','LineWidth',3,'Color',[0 0 0]);
ylabel('p'' (Pa)')
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
xlabel('time (hours)'); 

% quarter-phase markers. 90 deg of phase is (pi/2)/om seconds, so the spacing
% comes from the frequency and not from the hour axis -- change om and the
% lines follow. nq counts how many quarter cycles the record holds.
yyaxis left                             % so the labels sit on the velocity axis
nq = floor(om*t(end)/(pi/2) + 1e-9);
tq = (0:nq)*(pi/2)/om;                  % times of 0, 90, 180, ... deg [s]
for k = 1:numel(tq)
    xline(tq(k)/3600, ':', sprintf('%g\\circ', 90*(k-1)), ...
          'Color',[.5 .5 .5], 'LineWidth',0.75, 'FontSize',8, ...
          'LabelOrientation','horizontal', ...
          'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom');
end
xlim([0 t(end)/3600])
legend([h1 h2 h3],{'u (east)','v (north)','p'''},'Location','ne')
title('what a mooring would see')

%% --- two waves of equal amplitude from different directions -------------
% Nash et al. (2004, JPO 34, 1117-1134), Fig. 17: a wave with lambda_x = 80 km
% reflecting off a north-south wall makes a standing pattern in which the flux
% ACROSS the wall cancels everywhere, while the flux ALONG it reverses sign
% every quarter wavelength. A mooring sees none of that spatial structure --
% it sees one u, v, p record whose apparent flux direction depends on where in
% the pattern it happens to sit. Set A2 = A1+180 to reproduce that figure; any
% other pair of angles is the general two-wave case.
%
% Conventions are the cell above and synth_baroclinic's: for wave j,
%   theta_j = k_j . x - om*t + psi_j
%   u_par   =  U0*sin(theta_j)              along  e_j   = [sind Aj, cosd Aj]
%   u_perp  = -(f/om)*U0*cos(theta_j)       along  ep_j  = [-cosd Aj, sind Aj]
%   p'      =  P0*sin(theta_j)              in phase with u_par
% ep_j is 90 deg COUNTERclockwise from e_j, which with the minus sign on
% u_perp turns the vector clockwise -- the northern-hemisphere sense, same as
% the single-wave cell. Requires f, om, rho0, c1, U0 from that cell.

A1   = 215;                 % wave 1 propagation azimuth (toward), deg CW from N
A2   = 35;                  % wave 2 propagation azimuth (toward), deg CW from N
                            %   A2 = A1+180 -> Nash Fig. 17 (pure standing wave)
U0b  = U0;                  % SAME along-wave amplitude for both waves [m/s]
psi1 = 0;                   % phase of wave 1 at the mooring [deg]
psi2 = 90;                  % phase of wave 2 at the mooring [deg]
                            %   psi1-psi2 is the only thing that matters; it is
                            %   equivalently WHERE the mooring sits in the
                            %   interference pattern drawn in panel (d)

kmag = sqrt(om^2 - f^2)/c1;             % horizontal wavenumber [rad/m]
lam  = 2*pi/kmag;                       % wavelength [m]
P0b  = rho0*c1*U0b*sqrt(1 - (f/om)^2);  % pressure amplitude of EACH wave [Pa]

e1  = [ sind(A1)  cosd(A1)];   ep1 = [-cosd(A1)  sind(A1)];   % [east north]
e2  = [ sind(A2)  cosd(A2)];   ep2 = [-cosd(A2)  sind(A2)];
k1  = kmag*e1;                 k2  = kmag*e2;                 % [rad/m]

% ---- what the mooring records: the SUM of the two waves -----------------
tt  = linspace(0, 2*2*pi/om, 500);      % two M2 periods [s], as in cell 1
[u1,v1,p1] = wave_uvp(tt, om, f, U0b, P0b, e1, ep1, psi1);
[u2,v2,p2] = wave_uvp(tt, om, f, U0b, P0b, e2, ep2, psi2);
us = u1 + u2;   vs = v1 + v2;   ps = p1 + p2;

% ---- time-mean flux. Average over exactly one period (endpoint excluded,
% so the mean is unbiased) rather than over tt, whose endpoints double up.
ta = (0:999)/1000 * (2*pi/om);
[ua1,va1,pa1] = wave_uvp(ta, om, f, U0b, P0b, e1, ep1, psi1);
[ua2,va2,pa2] = wave_uvp(ta, om, f, U0b, P0b, e2, ep2, psi2);
Fnum = [mean((ua1+ua2).*(pa1+pa2))  mean((va1+va2).*(pa1+pa2))];   % [W/m^2]

% ---- the same flux in closed form. With phi_j = k_j.x + psi_j and
% D = phi_1 - phi_2, the self terms give (1/2)P0*U0*(e1+e2) -- the two waves
% carrying their own flux -- and the cross terms add
%   (1/2)P0*U0*[ cos(D)*(e1+e2) + (f/om)*sin(D)*(ep1-ep2) ].
% The cos(D) part rescales the vector sum; the sin(D) part points SIDEWAYS to
% both waves and is what Nash's black arrows are. For A2 = A1+180 the first
% part vanishes identically and only the sideways part survives -- and because
% it carries f/om, that residual is weak here (f/om = 0.14 at 7.5 deg N) where
% at Nash's 36.5 deg N it is 0.6.
Dm   = deg2rad(psi1 - psi2);            % D at the mooring (x = 0)
Fmoor = 0.5*P0b*U0b*( (1+cos(Dm))*(e1+e2) + (f/om)*sin(Dm)*(ep1-ep2) );
F1    = 0.5*P0b*U0b*e1;                 % wave 1 alone
F2    = 0.5*P0b*U0b*e2;                 % wave 2 alone
Fsum  = F1 + F2;                        % vector sum, i.e. no interference

az   = @(F) mod(atan2d(F(1),F(2)), 360);            % [east north] -> deg CW of N
mag  = @(F) hypot(F(1),F(2));

fprintf('\n--- two waves, |U| = %.3f m/s each, lambda = %.1f km ---\n', U0b, lam/1e3);
fprintf('wave 1 alone   : %.3f W/m^2 toward %6.1f deg\n', mag(F1),   az(F1));
fprintf('wave 2 alone   : %.3f W/m^2 toward %6.1f deg\n', mag(F2),   az(F2));
fprintf('vector sum     : %.3f W/m^2 toward %6.1f deg   (no interference)\n', ...
        mag(Fsum), az(Fsum));
if mag(Fmoor) < 1e-9*P0b*U0b
fprintf('AT THE MOORING : no net flux (standing wave, this phase)\n');
else
fprintf('AT THE MOORING : %.3f W/m^2 toward %6.1f deg   (what you would infer)\n', ...
        mag(Fmoor), az(Fmoor));
end
fprintf('closed form vs time average: %.2e W/m^2 difference\n', norm(Fmoor-Fnum));

% ---- interference pattern. D varies only along k1-k2, so the bands run
% perpendicular to that and repeat every Li.
dk = k1 - k2;
if norm(dk) < 1e-12*kmag                % identical directions: no pattern
    Li = lam;   uniform = true;
else
    Li = 2*pi/norm(dk);   uniform = false;
end
fprintf('interference bands every %.1f km', Li/1e3);
if uniform, fprintf('  (none -- the two waves are collinear)'); end
fprintf('\n');

figure('Position',[100 100 1000 760]);

% (a) the two propagation directions and the flux the mooring would report
subplot(2,2,1)
q1 = quiver(0,0, e1(1), e1(2), 0, 'LineWidth',2, 'MaxHeadSize',0.8, ...
            'Color',[0 0.447 0.741]); hold on
q2 = quiver(0,0, e2(1), e2(2), 0, 'LineWidth',2, 'MaxHeadSize',0.8, ...
            'Color',[0.850 0.325 0.098]);
fh = Fmoor/max(mag(Fmoor), eps);        % unit vector, only the angle matters
q3 = quiver(0,0, fh(1), fh(2), 0, 'k', 'LineWidth',3, 'MaxHeadSize',0.8);
sh = Fsum/max(mag(Fsum), eps);
q4 = quiver(0,0, sh(1), sh(2), 0, '--', 'LineWidth',1.5, 'MaxHeadSize',0.8, ...
            'Color',[.5 .5 .5]);
axis equal; axis([-1.2 1.2 -1.2 1.2]); grid on
xlabel('lon'); ylabel('lat')
legend([q1 q2 q3 q4], {sprintf('wave 1 -> %g\\circ',A1), ...
                       sprintf('wave 2 -> %g\\circ',A2), ...
                       sprintf('net flux -> %.0f\\circ', az(Fmoor)), ...
                       'vector sum'}, 'Location','southoutside')
title('direction each wave goes, and the flux you would infer')

% (b) what the mooring records -- same layout as the single-wave cell
subplot(2,2,2)
yyaxis left
h1 = plot(tt/3600, us, '-', 'LineWidth',1.5, 'Color',[0.850 0.325 0.098]); hold on
h2 = plot(tt/3600, vs, '-', 'LineWidth',1.5, 'Color',[0 0.447 0.741]);
ylabel('velocity (m s^{-1})')
yyaxis right
h3 = plot(tt/3600, ps, '-', 'LineWidth',3, 'Color',[0 0 0]);
ylabel('p'' (Pa)')
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
xlabel('time (hours)')
yyaxis left
nq2 = floor(om*tt(end)/(pi/2) + 1e-9);
tq2 = (0:nq2)*(pi/2)/om;
for k = 1:numel(tq2)
    xline(tq2(k)/3600, ':', sprintf('%g\\circ', 90*(k-1)), ...
          'Color',[.5 .5 .5], 'LineWidth',0.75, 'FontSize',8, ...
          'LabelOrientation','horizontal', ...
          'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom');
end
xlim([0 tt(end)/3600])
legend([h1 h2 h3], {'u (east)','v (north)','p'''}, 'Location','ne')
title('what the mooring sees (sum of the two)')

% (c) hodograph. One clean wave traces an f/om ellipse that turns clockwise;
% two of them interfering flatten it, and the sense can even reverse.
subplot(2,2,3)
one = tt <= 2*pi/om;
plot(u1(one), v1(one), '-',  'LineWidth',1.2, 'Color',[0 0.447 0.741]); hold on
plot(u2(one), v2(one), '-',  'LineWidth',1.2, 'Color',[0.850 0.325 0.098]);
plot(us(one), vs(one), 'k-', 'LineWidth',2);
plot(us(1), vs(1), 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);
axis equal; grid on
xlabel('u (east, m s^{-1})'); ylabel('v (north, m s^{-1})')
legend({'wave 1','wave 2','sum','t = 0'}, 'Location','best')
title('current ellipse over one M_2 period')

% (d) Nash Fig. 17. Shading is the normalized pressure variance <p'^2>/P0^2 =
% 1+cos(D) (the HKE/APE bands); arrows are the time-mean flux from the closed
% form above, evaluated on the grid. The mooring (red circle) samples one point
% of this, which is Nash's warning: u'p' has to be averaged over at least half
% an interference wavelength before it means what you think it means.
subplot(2,2,4)
xg = linspace(-Li, Li, 161);   yg = xg;
[Xg, Yg] = meshgrid(xg, yg);
Dg  = dk(1)*Xg + dk(2)*Yg + Dm;
Eg  = 1 + cos(Dg);
Fug = 0.5*P0b*U0b*( (e1(1)+e2(1))*(1+cos(Dg)) + (f/om)*(ep1(1)-ep2(1))*sin(Dg) );
Fvg = 0.5*P0b*U0b*( (e1(2)+e2(2))*(1+cos(Dg)) + (f/om)*(ep1(2)-ep2(2))*sin(Dg) );
pcolor(xg/1e3, yg/1e3, Eg); shading flat; hold on
colormap(gca, [linspace(1,0.70,64)' linspace(1,0.79,64)' linspace(1,0.92,64)']);
clim([0 2]);                            % Eg = 1+cos(D) spans [0 2] by construction
cb = colorbar; cb.Label.String = '\langle p''^2\rangle / P_0^2';
is = 1:12:numel(xg);
quiver(Xg(is,is)/1e3, Yg(is,is)/1e3, Fug(is,is), Fvg(is,is), 1.3, 'k', ...
       'LineWidth',1);
plot(0, 0, 'o', 'MarkerSize',9, 'LineWidth',2, 'MarkerEdgeColor',[0.6 0 0], ...
     'MarkerFaceColor','w');
axis equal; axis([xg(1) xg(end) yg(1) yg(end)]/1e3)
xlabel('east (km)'); ylabel('north (km)')
title(sprintf('mean flux through the pattern (bands every %.0f km)', Li/1e3))

function [u,v,p] = wave_uvp(t, om, f, U0, P0, ehat, eperp, psi_deg)
% One plane wave at a fixed point (s = 0), in the synth_baroclinic convention:
% theta = -om*t + psi, u_par = U0*sin(theta), u_perp = -(f/om)*U0*cos(theta).
th   = -om*t + deg2rad(psi_deg);
upar =  U0*sin(th);
uprp = -(f/om)*U0*cos(th);
u    = upar*ehat(1) + uprp*eperp(1);
v    = upar*ehat(2) + uprp*eperp(2);
p    = P0*sin(th);
end
