clc
clear
close all
% figure 

% Questa strategia parte dalla strategia standard e implementa il trasferimento 
% biellittico prima del cambio di piano per ottimizzare i costi di questa manovra 

% Tempo totale del trasferimento : 117.7117 h
% Costo propulsivo (DeltaV) : 4.7001 km/s
% Rapporto propulsivo - massa iniziale : 0.7975

% Costante gravitazionale
mu = 398600;

%% Parametri orbitali del punto sull'orbita iniziale

% Semiasse maggiore
a_i =  24400.0;
% Eccentricità
e_i = 0.728300;
% Inclinazione
i_i =  0.104700;
% RAAN - Ascensione retta del nodo ascendente
OM_i = 1.208000;
% Anomalia del pericentro
om_i = 3.107000;
% Anomalia vera 
theta_i =  2.462000;

[r_i, v_i] = par2car(a_i, e_i, i_i, OM_i, om_i, theta_i, mu);

h_i = cross(r_i, v_i);
h_i_norm = norm(h_i);
e_h_i = h_i/h_i_norm;

th0 = 0;
thf = 2*pi;
delta_th = 0.005;


%% Parametri orbitali del punto sull'orbita finale

% Posizione
r_x_f = 1894.036400;
r_y_f = -7182.528900;
r_z_f = 4175.118600;
r_f = [r_x_f, r_y_f, r_z_f];
% Velocità
v_x_f = 8.651000;
v_y_f = 0.544800;
v_z_f = -2.987000;
v_f = [v_x_f, v_y_f, v_z_f];

h_f = cross(r_f, v_f);
h_f_norm = norm(h_f);
e_h_f = h_f/h_f_norm;

[a_f, e_f, i_f, OM_f, om_f, theta_f] = car2par(r_f, v_f, mu);

th0 = 0;
thf = 2*pi;
delta_th = 0.005;

%% Cambio di piano
[~, ~, theta_plane, ~] = ...
    changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);

% Per valutare se sia più vantaggiosa la sola manovra per il cambio di piano
% oppure se lo sia quella a 3 impulsi con una strategia biellittica serve
% considerare l'andamento del costo in funzione della dimensione dell'ellisse
% ausiliaria

if(theta_plane > pi)
    theta_t = theta_plane - pi;
else
    theta_t = theta_plane + pi;
end

ra_f = a_f*(1 + e_f);
rp_f = a_f*(1 - e_f);

ra_i = a_i*(1 + e_i);
rp_i = a_i*(1 - e_i);

n = 10000;
ra_t = linspace(2*rp_i, 3*ra_f, n);
Delta_V_vect_plane_tot = zeros(n,1);
Delta_V_vect_1 = zeros(n,1);
Delta_V_vect_2 = zeros(n,1);
Delta_V_vect_circ_1 = zeros(n,1);
Delta_V_vect_circ_2 = zeros(n,1);
Delta_V_vect_plane = zeros(n,1);
Delta_V_vect_ecc = zeros(n,1);
DeltaT_1_vect = zeros(n,1);
DeltaT_2_vect = zeros(n,1);
DeltaT_vect = zeros(n,1);
a_1_vect = zeros(n,1);
e_1_vect = zeros(n,1);
om_1_vect = zeros(n,1);
om_2_vect = zeros(n,1);
a_2_vect = zeros(n,1);
e_2_vect = zeros(n,1);
theta_plane1_vect = zeros(n,1);
gamma_vect = zeros(n,1);
delta_vect = zeros(n,1);
DeltaV_tot_plane = 0;

for(j=1:n)
    [DeltaV_tot, DeltaV1, DeltaV2, DeltaV_circ_1, DeltaV_circ_2, ...
        DeltaV_plane, DeltaV_ecc, Deltat1, Deltat2, DeltaT, a_1, e_1, om_1, ...
        om_2, a_2, e_2, theta_plane1, gamma, delta, v] = ...
        biellipticTransferChangePlane(a_i, e_i, om_i, ...
        a_f, e_f, om_f, ra_t(j), mu, OM_i, i_i, OM_f, i_f, theta_plane);

    Delta_V_vect_plane_tot(j) = DeltaV_tot;
    Delta_V_vect_1(j) = DeltaV1;
    Delta_V_vect_2(j) = DeltaV2;
    Delta_V_vect_circ_1(j) = DeltaV_circ_1;
    Delta_V_vect_circ_2(j) = DeltaV_circ_2;
    Delta_V_vect_plane(j) = DeltaV_plane;
    Delta_V_vect_ecc(j) = DeltaV_ecc;
    DeltaT_1_vect(j) = Deltat1;
    DeltaT_2_vect(j) = Deltat2;
    DeltaT_vect(j) = DeltaT;
    a_1_vect(j) = a_1;
    e_1_vect(j) = e_1;
    a_2_vect(j) = a_2;
    e_2_vect(j) = e_2;
    om_1_vect(j) = om_1;
    om_2_vect(j) = om_2;
    theta_plane1_vect(j) = theta_plane1;
    gamma_vect(j) = gamma; % anomalia vera dell'orbita circolare in cui parte 
    % la manovra biellittica e quindi coincide con l'anomalia di pericentro 
    % della prima orbita ellittica
    delta_vect(j) = delta;

    if(j == 1 || DeltaV_tot_plane > Delta_V_vect_plane_tot(j))
        DeltaV_tot_plane = Delta_V_vect_plane_tot(j);
        k = j;
    end

end

% Disegno del pianeta
scatter3(r_i(1), r_i(2), r_i(3), 25, 'g', 'filled');
hold on

% Disegno dell'orbita iniziale e finale
plotorbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, delta_th, mu)
plotorbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, delta_th, mu)

if(v == 1)
    % Animazione dei trasferimenti
    plotorbit_improved(a_i, e_i, i_i, OM_i, om_i, theta_i, pi, delta_th, mu)
    plotorbit_improved(ra_i, 0, i_i, OM_i, 0, (om_i + pi), gamma, delta_th, mu)
    plotorbit_improved(a_1_vect(k), e_1_vect(k), i_i, OM_i, om_1_vect(k), 0, pi, delta_th, mu)
    plotorbit_improved(a_2_vect(k), e_2_vect(k), i_f, OM_f, om_2_vect(k), pi, 2*pi, delta_th, mu)
    plotorbit_improved(ra_f, 0, i_f, OM_f, 0, om_2_vect(k), (pi + om_f), delta_th, mu)
    plotorbit_improved(a_f, e_f, i_f, OM_f, om_f, pi, theta_f, delta_th, mu)
end

if(v == 2)
    % Animazione dei trasferimenti
    plotorbit_improved(a_i, e_i, i_i, OM_i, om_i, theta_i, pi, delta_th, mu)
    plotorbit_improved(ra_i, 0, i_i, OM_i, 0, (om_i + pi), gamma, delta_th, mu)
    plotorbit_improved(a_1_vect(k), e_1_vect(k), i_i, OM_i, om_1_vect(k), 0, pi, delta_th, mu)
    plotorbit_improved(a_2_vect(k), e_2_vect(k), i_i, OM_i, om_2_vect(k), pi, 2*pi, delta_th, mu)
    plotorbit_improved(ra_f, 0, i_f, OM_f, 0, om_2_vect(k), (pi + om_f), delta_th, mu)
    plotorbit_improved(a_f, e_f, i_f, OM_f, om_f, pi, theta_f, delta_th, mu)
end

if(v == 3)
    % Animazione dei trasferimenti
    plotorbit_improved(a_i, e_i, i_i, OM_i, om_i, theta_i, pi, delta_th, mu)
    plotorbit_improved(ra_i, 0, i_i, OM_i, 0, (om_i + pi), gamma, delta_th, mu)
    plotorbit_improved(a_1_vect(k), e_1_vect(k), i_f, OM_f, om_1_vect(k), 0, pi, delta_th, mu)
    plotorbit_improved(a_2_vect(k), e_2_vect(k), i_f, OM_f, om_2_vect(k), pi, 2*pi, delta_th, mu)
    plotorbit_improved(ra_f, 0, i_f, OM_f, 0, om_2_vect(k), (pi + om_f), delta_th, mu)
    plotorbit_improved(a_f, e_f, i_f, OM_f, om_f, pi, theta_f, delta_th, mu)
end

legend('Start of scenario','Initial orbit','Parking orbit', ...
    'Approach to first maneuver','First circularization','First circular orbit', ...
    'Start of bi-elliptic maneuver','First semi-ellipse of bi-elliptic transfer', ...
    'Plane change and second bi-elliptic maneuver','Second semi-ellipse of bi-elliptic transfer', ...
    'Second circularization','Second circular orbit','Ellipticization','End of scenario')

% === Plot Terra ===
r_terra = 6378;
[xe, ye, ze] = sphere(1000);
xe = xe * r_terra;
ye = ye * r_terra;
ze = ze * r_terra;
earthImg = imread('earth.jpg');
earthImg = flipud(earthImg);
earthImg = fliplr(earthImg);
earthImg = imresize(earthImg, size(xe));

h_earth = surf(xe, ye, ze, earthImg, ...
    'FaceColor', 'texturemap', ...
    'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', ...
    'SpecularStrength', 0.5);

% Rimuove la Terra dalla legenda
set(get(get(h_earth,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Calcolo dei tempi per lo scenario
delta_t_i = TOF(a_i, e_i, theta_i, pi, mu);
delta_t_f = TOF(a_f, e_f, pi, theta_f, mu);
delta_t_TOT = delta_t_i + DeltaT_vect(k) + delta_t_f;
flight_hours = delta_t_TOT/3600
flight_days = delta_t_TOT/(3600*24)

% Costo complessivo della missione
DeltaV_TOT = DeltaV_tot_plane 

figure
plot(ra_t(1:end), Delta_V_vect_plane_tot(1:end),'LineWidth', 1.5)
hold on
scatter(ra_t(k), Delta_V_vect_plane_tot(k))
text(ra_t(k), Delta_V_vect_plane_tot(k), '  \leftarrow Optimum', 'Color', 'red');
xlabel('Apogee of the bi-elliptic transfer [km]')
ylabel('Total transfer cost [km/s]')
title('Total \DeltaV cost as a function of r_{a,t}')
legend('Total \DeltaV')

figure
plot(ra_t(1:end), Delta_V_vect_plane(1:end),'LineWidth', 1.5)
hold on
scatter(ra_t(k), Delta_V_vect_plane(k))
plot(ra_t(1:end), Delta_V_vect_1(1:end),'LineWidth', 1.5)
scatter(ra_t(k), Delta_V_vect_1(k))
plot(ra_t(1:end), Delta_V_vect_2(1:end),'LineWidth', 1.5)
scatter(ra_t(k), Delta_V_vect_2(k))
plot(ra_t(1:end), Delta_V_vect_circ_2(1:end),'LineWidth', 1.5)
scatter(ra_t(k), Delta_V_vect_circ_2(k))
xline(ra_t(k), '--r', 'ra_t = 7.5413e+04 km', 'LineWidth', 1);
xlabel('Apogee of the bi-elliptic transfer [km]')
ylabel('Transfer cost [km/s]')
title('Trend of \DeltaV cost as a function of r\_a\_t')

legend('\DeltaV plane change','Optimal plane change','\DeltaV first bi-elliptic maneuver', ...
    'Optimal first bi-elliptic maneuver','\DeltaV second bi-elliptic maneuver', ...
    'Optimal second bi-elliptic maneuver','\DeltaV second circularization', ...
    'Optimal second circularization')

% DeltaV1_opt = 1.476747806430367e-04 km/s
% DeltaV2_opt = 0.351875636603173 km/s
% DeltaV_circ_1_opt = 1.263092537683692 km/s
% DeltaV_circ_2_opt = 0.407720684669953 km/s
% DeltaV_plane_opt = 1.205445231956458 km/s
% DeltaV_ecc_opt = 1.471886478420781 km/s

% ra_t_opt = 7.541365418236889e+04 km;
% Apogeo più grande rispetto a quello dell'orbita finale

% Dunque il costo della strategia è particolarmente condizionato dalla
% prima manovra di circolarizzazione, dall'ultima ellitticizzazione e dal
% cambio di piano (tutte cose necessarie anche nella strategia alternativa
% di circolarizzazione) 
