clc
clear
% close all
figure

% Strategia standard:
%
% 1. Cambio di piano
% 2. Cambio di pericentro
% 3. Cambio di forma

% Tempo totale del trasferimento : - h
% Costo propulsivo (DeltaV) : 5.9966 km/s
% Rapporto propulsivo - massa iniziale : 0.8697

% Strategia standard peggiore

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
[deltaV_plane, om_t, theta_plane, alpha] = ...
    changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f, mu);

% plotorbit(a_i, e_i, i_f, OM_f, om_t, th0, thf, delta_th, mu)

%% Cambio di anomalia di pericentro
[deltaV_peri, theta_1, theta_2] = changePericenterArg(a_i, e_i, om_t, om_f, mu);

% plotorbit(a_i, e_i, i_f, OM_f, om_f, th0, thf, delta_th, mu)

%% Cambio di forma
% Si valutano diverse strategie di trasferimeto per ricavare l'opzione più
% vantaggiosa per i costi della missione

% Approccio con trasferimento ellittico bitangente 
type = ['pa'; 'ap'; 'aa'; 'pp'];
delta_v_tot = zeros(2,1);
delta_T = zeros(2,1);
a_trasf_vect = zeros(2,1);
e_trasf_vect = zeros(2,1);

for(j=1:2)
    [delta_v_1, delta_v_2, delta_t, a_trasf, e_trasf] = ...
        bitangentTransfer (a_i, e_i, a_f, e_f, type(j,1:2), mu);
    delta_v_tot(j,1) = delta_v_1 + delta_v_2;
    delta_T(j,1) = delta_t;
    a_trasf_vect(j,1) = a_trasf;
    e_trasf_vect(j,1) = e_trasf;
end

deltaV_tot_bitan = min(delta_v_tot);

% Approccio con trasferimento biellittico
% ra_f = a_f*(1 + e_f);
% rp_i = a_i*(1 - e_i);
% 
% n = 10000;
% ra_t = linspace(ra_f + rp_i, 2*ra_f + rp_i, n);
% Delta_V_tot = zeros(n,1);
% Delta_T_biel = zeros(n,1);
% 
% for(j=1:n)
%     [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2, DeltaT] = ...
%         biellipticTransfer(a_i, e_i, a_f, e_f, ra_t(j), mu);
%     Delta_V_tot(j) = DeltaV1 + DeltaV2 + DeltaV3;
%     Delta_T_biel(j) = DeltaT;
% end
% 
% Delta_V_tot_biell = min(Delta_V_tot);
% 
% plot(ra_t,Delta_V_tot)
% Dall'analisi si evince che il trasferimento bitangente è molto più
% vantaggioso a livello di costo della manovra

%% Dunque si plotta l'ultimo trasferimento per raggiungere l'orbita finale

% Costo complessivo della missione
DeltaV_TOT = deltaV_plane + deltaV_peri + deltaV_tot_bitan

%% Parte Grafica

scatter3(r_i(1), r_i(2), r_i(3), 20, 'g', 'filled');
hold on
scatter3(r_x_f, r_y_f, r_z_f, 20, 'b', 'filled');
hold on

% Disegno del pianeta
scatter3(0, 0, 0, 100, 'b', 'filled');
hold on

% Disegno dell'orbita iniziale e finale
plotorbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, delta_th, mu)
plotorbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, delta_th, mu)

% Animazione dei trasferimenti
plotorbit_improved(a_i, e_i, i_i, OM_i, om_i, theta_i, theta_plane, delta_th, mu)
plotorbit_improved(a_i, e_i, i_f, OM_f, om_t, theta_plane, theta_1(2), delta_th, mu)
plotorbit_improved(a_i, e_i, i_f, OM_f, om_f, theta_2(2), pi, delta_th, mu)
plotorbit_improved(a_trasf_vect(2), e_trasf_vect(2), i_f, OM_f, om_f, pi, 0, delta_th, mu)
plotorbit_improved(a_f, e_f, i_f, OM_f, om_f, 0, theta_f, delta_th, mu)

