clc
close all
clear
% figure

% Costante gravitazionale solare
mu_s = 1.327e11;
% Costante gravitazionale terrestre (SOI_T)
mu = 398600;

% Unità astronomica
UA = 149597870.7;

% Massa della Terra
m_T =  5.974e24;
% Massa del Sole
M_S = 1.989e30;

%% Parametri orbitali dell'orbita eliocentrica della Terra

% Semiasse maggiore
a_T = 1.4946e08;
% Eccentricità
e_T = 0.016;
% Inclinazione
i_T = 9.1920e-05;
% RAAN - Ascensione retta del nodo ascendente
OM_T = 2.7847;
% Anomalia del pericentro
om_T = 5.2643;
% Anomalia vera iniziale
th0_T = 0;
thf_T = 2*pi;

delta_th = 0.005;
delta_t = 10000;

% plotorbit(a_T, e_T, i_T, OM_T, om_T, th0, thf, delta_th, mu_s)

p_T = a_T*(1 - e_T^2); % semilato retto dell'orbita terrestre

%% Parametri orbitali dell'orbita eliocentrica di trasferimento del satellite

% Dati del trasferimento ottimale dello scenario 2 standard
a_t_opt = 2.024805471260838e+08;
e_t_opt = 0.254483121458009;
i_t_opt = 0.032938968620602;
om_t_opt = 2.982926357953945;
OM_t_opt = 1.102187598272199;

theta1_t_opt = 0.155893443939461;
theta2_t_opt = 2.755132348139930;

thetaf_A = 2.792526803190927;
thetai_T = 2.475194211919231;

% % Semiasse maggiore
% a_t_opt = 1.848219275542337e+08;
% % Eccentricità
% e_t_opt = 0.185344860301687;
% % Inclinazione
% i_t_opt = 0.028484859392634; 
% % Anomalia del pericentro
% om_t_opt = 2.855993321445267; 
% % Ascensione retta del nodo ascendente
% OM_t_opt = 1.231800840855091; 
% % Anomalia inizio trasferimento
% theta1_t_opt = 0.282372259975278; 
% % Anomalia fine trasferimento
% theta2_t_opt = 4.752483335252830; 
% 
% % Anomalia vera terrestre di inizio trasferimento
% thetai_T = 2.602413028495237;
% % Anomalia vera dell'asteroide di fine trasferimento 
% thetaf_A = 4.383695398283848;

% Parametri cartesiani iniziali del trasferimento
[ri_t_vect, vt_vect_1] = par2car(a_t_opt, e_t_opt, i_t_opt, OM_t_opt, om_t_opt, theta1_t_opt, mu_s);
[ri_T_vect, vT_vect_i] = par2car(a_T, e_T, i_T, OM_T, om_T, thetai_T, mu_s);

% Distanza della Terra dal Sole nel punto di inizio manovra
R_T_elio = p_T/(1 + e_T*cos(thetai_T));
% Calcolo della SOI terrestre nel punto di inizio trasferimento
r_SOI_T = R_T_elio*(m_T/M_S)^(2/5);

%% Parametri orbitali del punto sull'orbita di parcheggio (geocentrica)

% Posizione
r_x_i = 1894.036400;
r_y_i = -7182.528900;
r_z_i = 4175.118600;
r_i_vect = [r_x_i, r_y_i, r_z_i];
% Velocità
v_x_i = 8.651000;
v_y_i = 0.544800;
v_z_i = -2.987000;
v_i_vect = [v_x_i, v_y_i, v_z_i];

h_i_vect = cross(r_i_vect, v_i_vect);
h_i_norm = norm(h_i_vect);
e_h_i = h_i_vect/h_i_norm;

[a_i, e_i, ~, ~, ~, ~] = car2par(r_i_vect, v_i_vect, mu);

p_i = (1 - e_i^2)*a_i; % semilato retto

%% Calcolo dei parametri dell'orbita di fuga

% Calcolo dell'Orbita Iperbolica

% Velocità di fuga
v_h_inf = norm(vt_vect_1 - vT_vect_i);
% Velocità al perigeo dell'orbita di parcheggio
vp_i = sqrt(mu/p_i)*(1 + e_i);
% Pericentro dell'orbita ellittica (e iperbolica)
r_p = a_i*(1 - e_i);

% Data la velocità da raggiungere alla distanza r_SOI_T (-> infinito) e
% data la distanza al perigeo dell'orbita di fuga iperbolica, si può
% definire univocamente l'orbita attreverso i suoi parametri

% Si sfrutta la costante energetica per trovare il semiasse maggiore
% E = -mu/(2*a) = (1/2)*v_inf^2 - mu/r_SOI_T = (1/2)*v_inf^2

a = -mu/(v_h_inf^2); % Semiasse maggiore dell'iperbole (0 > a)
a_h = a;

% Dal semiasse maggiore è possibile ricavare l'eccentricità dell'orbita
% dalle relazioni che sussistono tra eccentricità-semiasse maggiore-perigeo

% a = (r_a + r_p)/2

% e = (r_a + r_p)/(r_a + r_p)
% r_a = r_p*((1+e)/(1-e))

% Sostituendo
% a = r_p/(1-e)

% e = 1 - r_p/a

e = 1 - r_p/a; % eccentricità dell'iperbole (e > 1)
e_h = e;

% Dunque si possono sfruttare le caratteristiche dell'iperbole per ricavare
% i restanti parametri (semiasse minore [b], asindoti, semidistanza focale
% [c], angolo di deflessione [delta_2] e parametro di impatto [DELTA])

% La semidistanza focale sarà:
c = r_p - a;
% Alternativamente e = c/a --> c = e*a
c_1 = -(e*a);

% Controllo
if(c_1 ~= c)
    fprintf('problema con i parametri orbitali')
end
              
%                        x^2   y^2
% Equazione iperbole ->  ___ - ___ = 1
%                        a^2   b^2

% L'asindoto si trova quando si fa tendere ad infinito l'ascissa della
% funzione:
% f(x) = y = +- sqrt(x^2*(b^2/a^2) - 1)

%  lim   +- sqrt(x^2*(b^2/a^2) - 1) =  lim  +- sqrt(x^2*(b^2/a^2))
% x->inf                       x->inf

% Asindoti -> y = +-(b/a)*x

% Dunque si può ricavare il semiasse minore:
b = sqrt(c^2 - a^2);

% L'angolo di deflessione è sin(delta_2) = a/c = 1/e
delta = 2*asin(1/e);

% Parametro di impatto (si trova come distanza tra gli asindoti e il
% fuoco, ovvero la Terra):

% Si considera il fascio di rette perpendicolari all'asindoto y = (b/a)*x:
% y = -(a/b)*x + q

% La retta ricercata deve passare dal fuoco -> F(c,0)
% q = (a/b) * c

% Dunque la retta cercata è : y = -(a/b)*(x - c)

% Si calcola il punto P(x_P,y_P) di intersezione tra le due rette:
% y = -(a/b)*(x - c) 
% y = (b/a)*x         
% Dunque si risolve per x:  -(a/b)*(x - c) = (b/a)*x
% Ottenendo : x_P = a^2*(e/c) --->  y_P = b/c*(a*e)

% Si calcola la distanza tra F e P:
% DELTA = dist(F,P) = sqrt((x_P - x_C)^2 + (y_P - y_C)^2)
% Questa distanza sarà il parametro di impatto:
DELTA = sqrt((c - a^2*(e/c))^2 + ((a*e)*(b/c))^2);

% Oppure
DELTA_1 = -a_h*sqrt(e_h^2 - 1);

% Infine calcoliamo la velocità richiesta al perigeo per arrivare alla
% velocità di fuga richiesta alla distanza limite della SOI terrestre)

% Ancora una volta sfruttiamo la conservazione dell'energia totale
% E = (1/2)*v^2 - mu/r = cost
% (1/2)*v_h_inf^2 - mu/r_SOI_T = (1/2)*vp_h^2 - mu/r_p
% Ma si ricorda che si fa tendere ad infinito il raggio della SOI terrestre
% nel sistema di riferimento geocentrico quindi: mu/r_SOI_T -> 0

vp_h = sqrt(v_h_inf + 2*mu/r_p); % velocità dell'orbita iperbolica al perigeo

deltaV = vp_h - vp_i; % Costo della manovra di fuga

% th0_h = 0;
% thf_h = acos((p_h/(2*a_i) - 1)/e_h);
% n = round((thf - th0)/delta_th);

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

% plotorbit(a_i, e_i, i_i, OM_i, om_i, theta_i, theta_i + 2*pi, delta_th, mu)

%% Parametri orbitali del punto sull'orbita finale (orbita di parcheggio)

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

[a_f, e_f, i_f, OM_f, om_f, theta_f] = car2par(r_f, v_f, mu);

th0 = 0;
thf = 2*pi;
delta_th = 0.005;

p_f = (1 - e_f^2)*a_f; % semilato retto

% plotorbit(a_f, e_f, i_f, OM_f, om_f, theta_f, theta_f + 2*pi, delta_th, mu)

%% RIPRESA DEL PRIMO SCENARIO
% Si inserisce anche il primo scenario per collegarlo con la manovra di
% fuga

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
scatter3(0, 0, 0, 100, 'blue', 'filled');
hold on

% Disegno del punto
scatter3(r_i(1), r_i(2), r_i(3), 25, 'g', 'filled');
hold on

p_h = a_h*(1 - e_h^2); 
th0_h_1 = 0;
thf_h_1 = acos((p_h/(3*a_f) - 1)/e_h);

% Disegno dell'orbita iniziale e finale
plotorbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, delta_th, mu)
plotorbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, delta_th, mu)
% Disegno dell'orbita di fuga
plotorbit_hyperbolic(a_h, e_h, i_f, OM_f, om_f, th0_h_1, thf_h_1, delta_th, mu)

% Animazione dei trasferimenti
plotorbit_improved(a_i, e_i, i_i, OM_i, om_i, theta_i, pi, delta_th, mu)
plotorbit_improved(ra_i, 0, i_i, OM_i, 0, (om_i + pi), theta_circ_1_ell, delta_th, mu)
plotorbit_improved(a_t, e_t, i_i, OM_i, theta_circ_1_ell, 0, pi, delta_th, mu)
plotorbit_improved(ra_f, 0, i_f, OM_f, 0, theta_circ_2, (om_f + pi), delta_th, mu)
plotorbit_improved(a_f, e_f, i_f, OM_f, om_f, pi, 2*pi, delta_th, mu)

th0_h = 0;
thf_h = acos((p_h/(2.5*a_f) - 1)/e_h);

th0_i = 0;
thf_i = 2*pi;

% Animazione della manovra di fuga
plotorbit_improved_hyperbolic(a_h, e_h, i_f, OM_f, om_f, th0_h, thf_h, delta_th, mu)

legend('Earth', 'Starting point', 'Initial orbit', 'Parking orbit', ...
    'Hyperbolic escape trajectory', 'Transfer from starting point to circularization start point', ...
    'First impulsive maneuver', 'First circular orbit', 'Second impulsive maneuver', ...
    'Connecting ellipse between circular orbits', 'Third impulsive maneuver', ...
    'Second circular orbit', 'Fourth impulsive maneuver', 'Transfer on final orbit', ...
    'Escape maneuver', 'Hyperbolic escape trajectory')


