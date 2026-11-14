clc
close all
clear
% figure

% Ottimizzazione della manovra di cattura del satellite nella SOI
% dell'asteroide

% deltaV_min = 0.234135152513084

% Costante gravitazionale solare
mu_s = 1.327e11;

m_A = 3.1464e10; % Massa dell'asteroide (kg)
G = 6.67e-11; % Costante di gravitazione universale
mu_A = G*m_A; % Costante gravitazionale asteroide (SOI_T)

% Unità astronomica
UA = 149597870.7;

% Massa del Sole (kg)
M_S = 1.989e30;

% Diametro dell'asteroide (km)
D = 0.33;

%% Parametri orbitali dell'orbita eliocentrica dell'asteroide

r_A_i_vect = [1.257283091380454e+08;
              2.179025686550188e+08;
             -1.109788941692904e+06];

v_A_i_vect = [-17.36883104272374;
              8.840929257352350;
              0.5486875533296964];

[a_A, e_A, i_A, OM_A, om_A, th0_A] = car2par(r_A_i_vect, v_A_i_vect, mu_s); % th0_A = 3.2753017708 rad
thf_A = th0_A + 2*pi;

% % Semiasse maggiore
% a_A = 1.324170*UA;
% % Eccentricità
% e_A = 0.280134;
% % Inclinazione
% i_A = 1.62*(pi/180);
% % RAAN - Ascensione retta del nodo ascendente
% OM_A = 69.08*(pi/180);
% % Anomalia del pericentro
% om_A = 162.84*(pi/180);

p_A = a_A*(1 - e_A^2);

delta_th = 0.005;

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

% Dati del trasferimento ottimale dello scenario 2 avanzato
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

% % Anomalia vera terrestre di inizio trasferimento
% thetai_T = 2.602413028495237;
% % Anomalia vera dell'asteroide di fine trasferimento 
% thetaf_A = 4.383695398283848;

%% Calcolo dei parametri utili per la fine del trasferimento del satellite

% Parametri cartesiani finali del trasferimento
[rf_t_vect, vt_vect_2] = par2car(a_t_opt, e_t_opt, i_t_opt, OM_t_opt, om_t_opt, theta2_t_opt, mu_s);
[rf_A_vect, vA_vect_f] = par2car(a_A, e_A, i_A, OM_A, om_A, thetaf_A, mu_s);

% Distanza dell'asteroide dal Sole nel punto di inizio manovra
R_A_elio = p_A/(1 + e_A*cos(thetaf_A));
% Calcolo della SOI dell'asteroide nel punto di fine trasferimento
r_SOI_A = R_A_elio*(m_A/M_S)^(2/5);

%% Ottimizzazione della cattura

% La cattura avviene nel punto di perigeo dell'orbita iperbolica (che
% diventa il punto di perigeo dell'orbita di cattura ellittica / circolare)

% Ma la velocità al perigeo dell'orbita di cattura condiziona il costo
% della manovra: tanto è più veloce l'orbita di cattura al perigeo
% (tanto più è eccentrica l'orbita e, quindi, vicina all'attrattore) tanto
% minore sarà il costo della manovra

% Sapendo la velocità d'ingresso nella SOI dell'asteroide, si fa variare il
% parametro della distanza di pericentro e si calcolano i costi della
% manovra

% Velocità in ingresso
v_h_inf = norm(vt_vect_2 - vA_vect_f);

a_h = - mu_A/(v_h_inf^2); % Semiasse maggiore dell'iperbole (0 > a)

rp_h_vect = linspace(D/2+0.1, r_SOI_A, 1e3);
N = length(rp_h_vect);

e_f_vect = linspace(0, 0.999, 1e2);
n = length(e_f_vect);

t = 0;
deltaV_min = 0;
deltaV_map = NaN(N,n);

for(i=1:N)
    r_p = rp_h_vect(i);

    a = a_h;
    e = 1 - r_p/a; % eccentricità dell'iperbole (e > 1)

    % La semidistanza focale sarà:
    c = r_p - a;

    b = sqrt(c^2 - a^2); % Semiasse minore

    % L'angolo di deflessione è sin(delta_2) = a/c = 1/e
    delta = 2*asin(1/e);

    % Parametro di impatto
    % DELTA = sqrt((c - a^2*(e/c))^2 + ((a*e)*(b/c))^2);
    DELTA = - a*sqrt(e^2 - 1);

    if(DELTA > D/2)
        % Velocità dell'orbita iperbolica al perigeo
        vp_h = sqrt(v_h_inf + 2*mu_A/r_p);

        % A questo punto si vanno a definire i parametri dell'ellisse di
        % cattura ottimale facendo variare il parametro di eccentricità (il
        % perigeo è già definito)

        for(j=1:n)   
            e_f = e_f_vect(j); % eccentricità ellisse di cattura
            
            % Si calcolano i parametri orbitali dell'ellisse di cattura
            a_f = r_p/(1 - e_f); % semiasse maggiore
            ra_f = r_p*(1 + e_f)/(1 - e_f); % distanza all'apogeo

            if(r_SOI_A > ra_f)
                p_f = a_f*(1 - e_f^2); % semilato retto
    
                % Velocità sull'orbita ellittica al perigeo
                vp_f = sqrt(mu_A/p_f)*(1 + e_f); 
    
                % Costo della manovra di fuga
                deltaV = vp_h - vp_f; 
                deltaV_map(i,j) = deltaV; 
        
                if(deltaV_min > deltaV || t == 0)
                    t = t+1;
                    deltaV_min = deltaV;
    
                    % Parametri ottimali dell'iperbole
                    e_h_opt = e;
                    rp_h_opt = r_p;
                    delta_h_opt = delta;
                    DELTA_h_opt = DELTA;
    
                    % Parametri ottimali dell'ellisse di cattura
                    e_f_opt = e_f;
                    a_f_opt = a_f;
                end
            end
        end
    end
end

% Plot dell'asteroide
hold on
scatter3(0, 0, 0, 50, 'b','filled')

p_h = a_h*(1 - e_h_opt^2); 
thf = 0;
th0_2 = - acos((p_h/(2*a_f_opt) - 1)/e_h_opt);
th0_1 = - acos((p_h/(2.5*a_f_opt) - 1)/e_h_opt);

% Plot delle orbite
plotorbit_hyperbolic(a_h, e_h_opt, 0, 0, 0, th0_1, thf, delta_th, mu_A)
plotorbit(a_f_opt, e_f_opt, 0, 0, 0, 0, 2*pi, delta_th, mu_A)

% Animazione dello scenario
plotorbit_improved_hyperbolic(a_h, e_h_opt, 0, 0, 0, th0_2, thf, delta_th, mu_A)
plotorbit_improved(a_f_opt, e_f_opt, 0, 0, 0, 0, 2*pi, delta_th, mu_A)

legend('Asteroid','Incoming hyperbolic orbit','Elliptic capture orbit', ...
    'Hyperbolic trajectory','Capture maneuver','Elliptic trajectory')

% Mashgrid dei parametri
[RP, EF] = meshgrid(rp_h_vect, e_f_vect);

figure
surf(RP/1000, EF, deltaV_map',  'EdgeColor', 'none')
xlabel('Hyperbola periapsis [km * 10e-3]')
ylabel('Capture orbit eccentricity [-]')
zlabel('\DeltaV [m/s]')
title('Trend of \DeltaV cost as a function of r_p and e_f')
colorbar;
view(135, 30) % Per vedere da una prospettiva gradevole

% Calcolo gradiente su tutta la griglia (una sola volta)
[dV_drp_map, dV_def_map] = gradient(deltaV_map', rp_h_vect, e_f_vect);

% Inizializzazione
rp = rp_h_vect(1);     % [km]
ef = e_f_vect(1);
alpha = 1.1e-2;    % lunghezza passo
max_iter = 70;
tol = 1e-6;
j = 1;
k = 1;

% Salvataggio percorso
path_rp = rp;
path_ef = ef;

for(i = 1:max_iter)  
    % Vettore gradiente nel punto 
    grad = [dV_drp_map(j,k); dV_def_map(j,k)];

    % Normalizzazione
    if norm(grad) == 0
        disp('Gradiente nullo: minimo trovato.');
        break;
    end
    % Direzione di massima crescita
    dir = grad/norm(grad);
    
    % Discesa (muovo nella direzione opposta al gradiente)
    rp_new = rp - alpha * dir(1);
    ef_new = ef - alpha * dir(2);
    % (direzione di massima decrescita)

    % Nuova posizone sulla mappa delle derivate parziali delle due
    % variabili
    j = j + round(-dir(2));
    k = k + round(-dir(1));
    
    % Controlla limiti dominio
    if rp_new < min(rp_h_vect) || rp_new > max(rp_h_vect) || ...
       ef_new < min(e_f_vect) || ef_new > max(e_f_vect)
        disp('Fuori dominio, interrotto');
        break;
    end

    % Verifica convergenza
    if norm([rp_new - rp, ef_new - ef]) < tol
        break;
    end
    
    % Aggiorna posizione
    rp = rp_new;
    ef = ef_new;
    
    % Salva percorso
    path_rp(end+1) = rp;
    path_ef(end+1) = ef;
end

% Plot finale
figure;
contourf(RP/1000, EF, deltaV_map', 40); 
hold on;
xlabel('Hyperbola periapsis [km]');
ylabel('Capture orbit eccentricity [-]');
title('Contours of \DeltaV cost');
colorbar;

% Percorso gradiente
plot(path_rp/1000, path_ef, 'r-', 'LineWidth', 2);
scatter(path_rp(1)/1000, path_ef(1), 50, 'k', 'filled'); ...
    text(path_rp(1)/1000, path_ef(1), '  Start');
scatter(rp_h_opt/1000, e_f_opt, 70, 'r', 'filled'); ...
    text(rp_h_opt/1000, e_f_opt, '  \leftarrow Optimum', 'Color', 'red');


% % Mashgrid dei parametri
% [RP, EF] = meshgrid(rp_h_vect, e_f_vect);
% 
% figure
% surf(RP/1000, EF, deltaV_map',  'EdgeColor', 'none')
% xlabel('Hyperbola periapsis [km * 10e-3]')
% ylabel('Capture orbit eccentricity [-]')
% zlabel('\DeltaV [m/s]')
% title('Trend of \DeltaV cost as a function of r_p and e_f')
% colorbar;
% view(135, 30) % Per vedere da una prospettiva gradevole
%
% figure;
% contourf(RP/1000, EF, deltaV_map', 40);  % 40 livelli di contorno
% xlabel('Hyperbola periapsis [km]');
% ylabel('Capture orbit eccentricity [-]');
% title('Contours of \DeltaV cost');
% colorbar;
% hold on;
% 
% % Evidenzia punto ottimo
% scatter(rp_h_opt/1000, e_f_opt, 100, 'r', 'filled');
% text(rp_h_opt/1000, e_f_opt, '  \leftarrow Optimum', 'Color', 'red');


