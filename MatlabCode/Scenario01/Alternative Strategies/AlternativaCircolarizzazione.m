clc
clear
close all
% figure

% Questa strategia parte dalla strategia standard e implementa una circolarizzazione
% prima del cambio di piano per ottimizzare i costi di questa manovra 

% Tempo totale del trasferimento : 35.2420 h
% Costo propulsivo (DeltaV) : 4.7001 km/s
% Rapporto propulsivo - massa iniziale : 0.7975

% Costante gravitazionale
mu = 398600;

%% Parametri orbitali del punto sull'orbita iniziale

% Semiasse maggiore
a_i = 24400.0;
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

% plotorbit(a_i, e_i, i_i, OM_i, om_i, theta_i, theta_i + 2*pi, delta_th, mu)

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

% plotorbit(a_f, e_f, i_f, OM_f, om_f, theta_f, theta_f + 2*pi, delta_th, mu)

%% Cambio di piano
[~, ~, theta_plane, ~] = ...
    changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);

% alpha_grad = alpha*(90/pi);
% 
% plotorbit(a_f, e_f, i_i, OM_i, om_t, th0, thf, delta_th, mu)

% Per valutare se sia più vantaggiosa la sola manovra per il cambio di piano
% oppure se lo sia quella a 3 impulsi con una strategia biellittica serve
% considerare l'andamento del costo in funzione della dimensione dell'ellisse
% ausiliaria

% if(theta_plane > pi)
%     theta_t = theta_plane - pi;
% else
%     theta_t = theta_plane + pi;
% end

if(0 >  theta_plane)
    theta_plane = theta_plane + 2*pi;
end

ra_f = a_f*(1 + e_f);
rp_f = a_f*(1 - e_f);

ra_i = a_i*(1 + e_i);
rp_i = a_i*(1 - e_i);

% ------------------------------------------------------------------
[deltaV_tot, deltaV_circ_1, deltaV_trasf, deltaV_circ_2, deltaV_ecc, deltaV_plane, ...
    deltaT_circ_1, deltaT_circ_2, deltaT_t, deltaT, a_t, e_t, om_t, theta_circ_1_ell, ...
    theta_circ_2, theta_circ_2_ell] = changePlaneCircTransfer(a_i, e_i, om_i, a_f, e_f, ...
    om_f, mu, OM_i, i_i, OM_f, i_f, theta_plane);
% ------------------------------------------------------------------

% Disegno del pianeta
scatter3(r_i(1), r_i(2), r_i(3), 25, 'b', 'filled');
hold on

% Disegno dell'orbita iniziale e finale
plotorbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, delta_th, mu)
plotorbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, delta_th, mu)

% Animazione dei trasferimenti
plotorbit_improved(a_i, e_i, i_i, OM_i, om_i, theta_i, pi, delta_th, mu)
plotorbit_improved(ra_i, 0, i_i, OM_i, 0, (om_i + pi), theta_circ_1_ell, delta_th, mu)
plotorbit_improved(a_t, e_t, i_i, OM_i, theta_circ_1_ell, 0, pi, delta_th, mu)
plotorbit_improved(ra_f, 0, i_f, OM_f, 0, theta_circ_2, (om_f + pi), delta_th, mu)
% plotorbit(ra_f, 0, i_f, OM_f, 0, (gamma + pi), gamma + (2*pi - (pi + om_f)), delta_th, mu)
plotorbit_improved(a_f, e_f, i_f, OM_f, om_f, pi, theta_f, delta_th, mu)

legend('Start of scenario','Initial orbit','Parking orbit', ...
    'Approach to first maneuver','First circularization','First circular orbit', ...
    'Start of bi-elliptic maneuver','First semi-ellipse of bi-elliptic transfer', ...
    'Plane change and circularization','Second circular orbit','Ellipticization', ...
    'Path on parking orbit')

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
deltaT_i = TOF(a_i, e_i, theta_i, pi, mu);
deltaT_f = TOF(a_f, e_f, pi, theta_f, mu);
deltaT_tot = deltaT + deltaT_i + deltaT_f;

flight_hours = deltaT_tot/3600
flight_days = deltaT_tot/(3600*24)

% Costo complessivo della missione
DeltaV_TOT = deltaV_tot 

% figure
% plot(ra_t, Delta_V_vect_plane)

% DeltaV_trasf = 0.407560260385891 km/s
% DeltaV_circ_1 = 1.471886478420781 km/s
% DeltaV_circ_2 = 0.351908094415124 km/s
% DeltaV_plane = 1.205699447516880 km/s
% DeltaV_ecc = 1.263092537683692 km/s

