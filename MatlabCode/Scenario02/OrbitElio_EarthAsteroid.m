clc
clear
close all
% figure

% Questo script ha l'unico scopo di mostrare l'animazione del moto dei tre
% corpi (Terra, asteroide e satellite) durante il trasferimento 
% interplanetario ottimo calcolato nello Scenario 02

% Costante gravitazionale solare
mu_s = 1.327e11;

% Unità astronomica
UA = 149597870.7;

%% Parametri orbitali dell'orbita eliocentrica della Terra

% Semiasse maggiore
a_T = 1.4946e08;
% Eccentricità
e_T = 0.016;
% Inclinazione
i_T =  9.1920e-05;
% RAAN - Ascensione retta del nodo ascendente
OM_T = 2.7847;
% Anomalia del pericentro
om_T = 5.2643;
% Anomalia vera iniziale
th0_T = 0;
thf_T = 2*pi;

delta_th = 0.005;
th0 = 0;
thf = 2*pi;
delta_t = 10000;

plotorbit(a_T, e_T, i_T, OM_T, om_T, th0, thf, delta_th, mu_s)

%% Parametri orbitali dell'orbita eliocentrica dell'asteroide

% Semiasse maggiore
a_A = 1.324170*UA;
% Eccentricità
e_A = 0.280134;
% Inclinazione
i_A = 1.62*(pi/180);
% RAAN - Ascensione retta del nodo ascendente
OM_A = 69.08*(pi/180);
% Anomalia del pericentro
om_A = 162.84*(pi/180);
% Anomalia vera iniziale
th0_A = 2*pi - 1.85;
thf_A = 1.85;

plotorbit(a_A, e_A, i_A, OM_A, om_A, th0, thf, delta_th, mu_s)

%% Plot del Sole
scatter3(0, 0, 0, 200, '*', 'o');
hold on

%% Plot del modo dei corpi celesti

a_t_opt = 1.958788230496195e+08;
e_t_opt = 0.228733651775526;
i_t_opt = 0.024165846383918; 
om_t_opt = 3.109859394462623; 
OM_t_opt = 0.957949503260746; 
theta1_t_opt = 0.028056658797491; 
theta2_t_opt = 4.306212641381219; 
deltaV_min = 4.870790424655530;
thetaf_A = 4.326337702283748;
thetai_T = 2.330051938232224;

plotcelestialbodies(a_T, e_T, i_T, OM_T, om_T, th0_T, thf_T, a_A, e_A, i_A, OM_A, ...
    om_A, th0_A, delta_t, mu_s, theta1_t_opt, theta2_t_opt, thetaf_A, thetai_T, ...
    a_t_opt, e_t_opt, i_t_opt, om_t_opt, OM_t_opt)

legend('Orbita terrestre','Orbita asteroide','Sole','Punto di partenza Terra',...
    'Punto di partenza asteroide','Traiettoria orbitale Terra','Traiettoria orbitale asteroide',...
    'Traiettoria di trasferimento satellite','Terra','Asteroide','Satellite')

N_tot = 4722; % Intervalli temporali di discretizzazione totali
N_sat = 2814; % Intervalli temporali di discretizzazione dalla partenza del satellite
% N_trasf = 2726; % Intervalli temporali di discretizzazione del trasferimento

deltaT_tot = 10000*N_tot; % Tempo totale animazione
deltaT_sat = 10000*N_sat; % Tempo dal lancio del satellite
% deltaT_trasf = 10000*N_trasf; % Tempo totale del trasferimento

% flight_hours = deltaT_trasf/3600
% flight_days = flight_hours/24
% flight_weeks = flight_days/7
% flight_years = flight_weeks/53

deltaV = deltaV_min

