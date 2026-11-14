clc
clear
close all

% Il secondo scenario si propone di ottimizzare un trasferimento diretto da
% un'orbita eliocentrica della Terra ad un'orbita eliocentrica di un
% asteroide assegnato 

% Questo script fa riferimento allo scenario 2 con rendezvous tra satellite 
% e asteroide con analisi temporale dal 4 gennaio 2025 al 3 gennaio 2028

% Tempo stimato per lo script: 6 ore

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
i_T = 9.1920e-05;
% RAAN - Ascensione retta del nodo ascendente
OM_T = 2.7847;
% Anomalia del pericentro
om_T = 5.2643;

th0_T = 0;
thf_T = 2*pi;

%% Parametri orbitali dell'orbita eliocentrica dell'asteroide

% % Semiasse maggiore
% a_A =  1.324170*UA;
% % Eccentricità
% e_A = 0.280134;
% % Inclinazione
% i_A =  1.62*(pi/180);
% % RAAN - Ascensione retta del nodo ascendente
% OM_A = 69.08*(pi/180);
% % Anomalia del pericentro
% om_A = 162.84*(pi/180);

% Calcolo dell'anomalia vera dell'asteroide corrispondente all'anomalia
% vera terrestre nulla (4 Gennaio 2025)

% I dati per questa parte sono stati presi dal sito 
% https://ssd.jpl.nasa.gov/horizons/app.html#/

% Dal sito è stato possibile individuare le coordinate del vettore
% posizione dell'asteroide nel sistema di riferimento eliocentrico e le 
% componenti del vettore velocità nel momento in cui la Terra orbita nel 
% suo punto di perigeo 

r_A_i_vect = [1.257283091380454e+08;
              2.179025686550188e+08;
             -1.109788941692904e+06];

v_A_i_vect = [-17.36883104272374;
              8.840929257352350;
              0.5486875533296964];

[a_A, e_A, i_A, OM_A, om_A, th0_A] = car2par(r_A_i_vect, v_A_i_vect, mu_s); % th0_A = 3.2753017708 rad
thf_A = th0_A + 2*pi;

p_A = a_A*(1 - e_A^2);

%% Ottimizzazione del trasferimento orbitale

tic;

% Periodo orbitale completo dell'asteroide 
T_A = TOF(a_A, e_A, 0, 2*pi, mu_s);
% Periodo orbitale terrestre
T_T = TOF(a_T, e_T, 0, 2*pi, mu_s);

% Passo di discretizzazione del tempo
dt = 10000;

% Vettore dei tempi
t_vect = 0 : dt : T_A;
t_vect_T = 0 : dt : T_T;

N = length(t_vect);
N_T = length(t_vect_T);

% Inizializzazione del vettore delle anomalie vere
theta_A_vect = zeros(1, N);
theta_A_vect(1) = th0_A;
theta_T_vect = zeros(1, N_T);
theta_T_vect(1) = th0_T;

% Inizializzazione dei vettori contenenti i vettori posizione e velocità
% per ottimizzare il ciclo for
r_vect_T_v = ones(3, N_T);
v_vect_T_v = ones(3, N_T);
r_vect_A_v = ones(3, N);
v_vect_A_v = ones(3, N);

% Il piano orbitale dell'orbita di trasferimento è dipendente dal punto di
% partenza sull'orbita terrestre e quello di arrivo sull'orbita
% dell'asteroide (theta1_i e theta2_f)

for(j = 1:N_T-1)
    % Calcolo dell'anomalia vera successiva
    [~, theta2] = KepleroDiretto(a_T, e_T, theta_T_vect(j), dt, mu_s);
    theta_T_vect(j+1) = theta2;

    % Calcolo del vettore corrispondente
    [r_vect_T, v_vect_T] = par2car(a_T, e_T, i_T, OM_T, om_T, theta2, mu_s);
    r_vect_T_v(:,j) = r_vect_T;
    v_vect_T_v(:,j) = v_vect_T;
end

for(j = 1:N-1)
    % Calcolo dell'anomalia vera successiva
    [~, theta2] = KepleroDiretto(a_A, e_A, theta_A_vect(j), dt, mu_s);
    theta_A_vect(j+1) = theta2;

    % Calcolo del vettore corrispondente
    [r_vect_A, v_vect_A] = par2car(a_A, e_A, i_A, OM_A, om_A, theta2, mu_s);
    r_vect_A_v(:,j) = r_vect_A;
    v_vect_A_v(:,j) = v_vect_A;
end

Time_1 = toc;
disp(['Tempo impiegato: ', num2str(Time_1), ' secondi']);

tic;

% Si prende come parametro di riferimento una finestra temporale
% opportunamente verificata per un'anomalia vera (in questo caso
% all'apogeo) e si assume che la finestra temporale associata ad un'altra
% anomalia vera sia proporzionale al rapporto tra la velocità
% di riferimento e la velocità associata alla suddetta anomalia

m_A = 3.1464e10; % Massa dell'asteroide (kg)
G = 6.67e-11; % Costante di gravitazione universale
mu_A = G*m_A; % Costante gravitazionale asteroide (SOI_T)

% Massa del Sole (kg)
M_S = 1.989e30;

% Raggio medio della SOI dell'asteroide Itokawa (serve per la finestra di
% tolleranza temporale dell'ottimizzazione)
R_A_elio_1 = p_A/(1 + e_A*cos(0));
r_SOI_1 = R_A_elio_1*(m_A/M_S)^(2/5);

R_A_elio_2 = p_A/(1 + e_A*cos(pi));
r_SOI_2 = R_A_elio_2*(m_A/M_S)^(2/5);

r_SOI_med = 0.5*(r_SOI_1 + r_SOI_2);

% Si calcola l'intervallo temporale di riferimento per il calcolo delle
% finestre temporali
Deltat_rif = 100*timeWindow(a_A, e_A, pi, r_SOI_med, mu_s);
% Si aumenta la finestra temporale di due grandezze di misura per rendere
% più congruente le differenze di discretizzazione temporale

Deltat = zeros(N,1); % Finestre temporali associate alle anomalie vere

p_A = a_A*(1 - e_A^2);
v = zeros(N,1);
v_rif = sqrt(mu_s/p_A)*(1 - e_A);

for(i = 1:N)
    v_r = sqrt(mu_s/p_A)*e_A*sin(theta_A_vect(i));
    v_t = sqrt(mu_s/p_A)*(1+e_A*cos(theta_A_vect(i)));
    v(i) = sqrt(v_r^2 + v_t^2);

    Deltat(i) = (v_rif/v(i))*Deltat_rif;
end

Time2 = toc;
disp(['Tempo impiegato: ', num2str(Time2), ' secondi']);

tic;

% Numero di anomalie di pericentro da considerare nell'ottimizzazione
n = 100;
% Anomalie di pericentro dell'orbita di trasferimento
om_t_vect = linspace(0, 2*pi, n); 

% Si procede con un'ottimizzazione dei costi con dei cicli for annidati che
% facciano variare i parametri delle anomalie vere (iniziale e finale) e
% dell'anomalia di pericentro dell'orbita di trasferimento

% Inizializzazione della variabile per l'ottimizzazione dei costi
deltaV_min = 0;
deltaV_map_1 = zeros(N_T,n);
deltaV_map_2 = zeros(N_T,n);
deltaV_map_3 = zeros(N_T,n);
deltaV_min_1 = 0;
deltaV_min_2 = 0;
deltaV_min_3 = 0;
deltaV_vect_1 = zeros(N_T,1);
deltaV_vect_2 = zeros(N_T,1);
deltaV_vect_3 = zeros(N_T,1);

% Inizializzazione per il conto delle rivoluzioni della Terra prima 
% dell'inizio del trasferimento ottimale
riv_T = 0;

% Ricerca del trasferimento ottimale prendendo come finestra temporale un
% intervallo pari a tre rivoluzioni della Terra
days = round(3*T_T/(3600*24)); % 1094 giorni considerati
% Inizio: 4 Gennaio 2025
% Fine: 3 Gennaio 2028

for(k = 1:3*N_T)

    % Selezione dei dati iniziali considerati
    if(k > N_T)
        if(k > 2*N_T)
            riv_T = 2;
            theta_T = theta_T_vect(k-2*N_T);
            r_vect_T = r_vect_T_v(:,k-2*N_T);
            v_vect_T = v_vect_T_v(:,k-2*N_T);
        else
            riv_T = 1;
            theta_T = theta_T_vect(k-N_T);
            r_vect_T = r_vect_T_v(:,k-N_T);
            v_vect_T = v_vect_T_v(:,k-N_T); 
        end
    else
        theta_T = theta_T_vect(k);
        r_vect_T = r_vect_T_v(:,k);
        v_vect_T = v_vect_T_v(:,k);
    end

    if(k > N)
        if(k > 2*N)
            % Anomalia vera iniziale dell'asteroide
            theta_A_i = theta_A_vect(k-2*N);
        else
            % Anomalia vera iniziale dell'asteroide
            theta_A_i = theta_A_vect(k-N);
        end
    else
        % Anomalia vera iniziale dell'asteroide
        theta_A_i = theta_A_vect(k);
    end

    for(s = 1:N)

        % Anomalia finale dell'asteroide
        if((s+k) > N)
            if((s+k) > 2*N)
                theta_A_f = theta_A_vect(s+k - 2*N);
                r_vect_A_f = r_vect_A_v(:, s+k - 2*N);
                v_vect_A_f = v_vect_A_v(:, s+k - 2*N);

                t = Deltat(s+k - 2*N);
            else
                theta_A_f = theta_A_vect(s+k - N);
                r_vect_A_f = r_vect_A_v(:, s+k - N);
                v_vect_A_f = v_vect_A_v(:, s+k - N);

                t = Deltat(s+k - N);
            end
        else
            theta_A_f = theta_A_vect(s+k);
            r_vect_A_f = r_vect_A_v(:, s+k);
            v_vect_A_f = v_vect_A_v(:, s+k);

            t = Deltat(s+k);
        end

        for(j = 1:n)

            % Anomalia di pericentro orbita di trasferimento
            om_t = om_t_vect(j);

            % Si possono calcolare tutti i parametri dell'orbita di
            % trasferimento considerando i costi della missione

            % Si definisce il piano orbitale di trasferimento che congiunge i
            % due punti dello spazio definiti dalle posizioni iniziale e finale
            % sulle due orbite (date dalle anomalie)
            vect = cross(r_vect_A_f, r_vect_T);
            h_t_vect = vect/norm(vect);
        
            % Versori direzionali della terna di riferimento
            i_vect = [1; 0; 0];
            j_vect = [0; 1; 0];
            k_vect = [0; 0; 1];
        
            % Inclinazione dell'orbita
            i_t = acos(dot(h_t_vect, k_vect));
        
            % Determinazione della linea dei nodi
            N_vect = cross(k_vect, h_t_vect)/norm(cross(k_vect, h_t_vect));
        
            % Calcolo dell'ascensione retta del nodo ascenendente 
            Ny = dot(N_vect,j_vect);
            Nx = dot(N_vect,i_vect);

            if(Ny>=0)
                OM_t = acos(Nx);
            else
                OM_t = 2*pi - acos(Nx);
            end
    
            % Matrici di rotazione
    
            R_om = [cos(OM_t)    sin(OM_t)     0;
                    -sin(OM_t)   cos(OM_t)     0;
                        0           0          1];
    
            R_i = [     1           0          0;
                        0        cos(i_t)   sin(i_t);
                        0       -sin(i_t)   cos(i_t)];
    
            R_w = [cos(om_t_vect(j)) sin(om_t_vect(j))           0;
                  -sin(om_t_vect(j)) cos(om_t_vect(j))           0;
                          0                  0                   1];
    
            T = R_w*R_i*R_om;
    
            % Definizione dei vettori posizione dei punti 1 e 2 sull'orbita 
            % di trasferimnto nel sistema perifocale
            rt_1_vect = T*r_vect_T;
            rt_2_vect = T*r_vect_A_f;
    
            % Definizione delle anomalie vere sull'orbita di trasferimento
            theta1_t_cos = acos(rt_1_vect(1)/norm(rt_1_vect));
            theta1_t_sin = asin(rt_1_vect(2)/norm(rt_1_vect));
    
            if(theta1_t_sin > 0)
                theta1_t = theta1_t_cos;
            else
                theta1_t = 2*pi - theta1_t_cos;
            end
            
            theta2_t_cos = acos(rt_2_vect(1)/norm(rt_2_vect));
            theta2_t_sin = asin(rt_2_vect(2)/norm(rt_2_vect));
    
            if(theta2_t_sin > 0)
                theta2_t = theta2_t_cos;
            else
                theta2_t = 2*pi - theta2_t_cos;
            end
    
            % Definizione dei parametri di forma della conica di
            % trasferimento tramite l'equazione della conica in coordinate
            % polari
    
            % Eccentricità dell'ellisse
            e_t = (norm(rt_1_vect) - norm(rt_2_vect))/(norm(rt_2_vect)*cos(theta2_t)...
                - norm(rt_1_vect)*cos(theta1_t));
    
            if((1 > e_t) && (e_t >= 0))
                % Semiasse maggiore
                a_t = norm(rt_1_vect)*((1 + e_t*cos(theta1_t))/(1 - e_t^2));
                % Raggio di pericentro 
                rp_t = a_t*(1 - e_t);

                % Si cerca di imporre una condizione per escludere tutte le
                % traiettorie orbitali che orbitano troppo vicino al Sole

                % Si impone una vicinanza minima di 0.7 UA, ma non si può
                % imporre attraverso la comparazione con il raggio di
                % perigeo perché questa condizione risulterebbe troppo
                % stringente per traiettorie di trasferimento che
                % coinvolgono solo una porzione di orbita vicino all'apogeo

                % Per ottimizzare i costi computazionali di questo calcolo
                % si fanno delle considerazioni preliminari per restringere
                % il campo di ricerca

                % Informazioni importanti sulla vicinante di un punto di
                % un'orbita da un attratore viene offerto dall'anomalia
                % vera

                if(theta2_t > theta1_t)
                    % In questo caso si considera la differenza di anomalia
                    % con 0 rad (per theta1_t) e con 2*pi rad (per
                    % theta2_t)

                    % Semilato retto
                    p_t = a_t*(1 - e_t^2);

                    if((2*pi - theta2_t) > theta1_t)
                        r_min = p_t/(1 + e_t*cos(theta1_t));
                    else
                        r_min = p_t/(1 + e_t*cos(theta2_t));
                    end
                else
                    % In questo caso l'orbita passa necessariamente per il
                    % perigeo
                    r_min = rp_t;
                end
           
                if(r_min > 0.7*UA)
                    % Calcolo dei tempi di trasferimento
                    delta_T_t = TOF(a_t, e_t, theta1_t, theta2_t, mu_s);
                    % Calcolo del tempo orbitale dell'asteroide dall'inizio del
                    % trasferimento
                    delta_T_A = TOF(a_A, e_A, theta_A_i, theta_A_f, mu_s);
    
                    if((delta_T_A + t > delta_T_t) && (delta_T_t > delta_T_A - t))

                        % Determinazione delle velocità nei punti della manovra
                        [~ , vt_vect_1] = par2car(a_t, e_t, i_t, OM_t, om_t, theta1_t, mu_s);
                        [~ , vt_vect_2] = par2car(a_t, e_t, i_t, OM_t, om_t, theta2_t, mu_s);
                
                        % Calcolo dei costi della manovra
                        deltaV1 = norm(vt_vect_1 - v_vect_T);
                        deltaV2 = norm(v_vect_A_f - vt_vect_2);

                        deltaV = deltaV1 + deltaV2;
                        if(k > N_T)
                            if(k > 2*N_T)
                                deltaV_map_3(k-2*N_T, j) = deltaV;
                                deltaV_vect_3(k-2*N_T) = deltaV;
                            else
                                deltaV_map_2(k-N_T, j) = deltaV; 
                                deltaV_vect_2(k-N_T) = deltaV;
                            end
                        else
                            deltaV_map_1(k, j) = deltaV;
                            deltaV_vect_1(k) = deltaV;
                        end
                
                        % Selezione del trasferimento ottimo
                        if(deltaV_min == 0 || deltaV_min > deltaV)
    
                            % Si applica una distinzione per capire il costo
                            % ottimale del traferimento nella specifica
                            % rivoluzione terrestre considerata
                            if(k > N_T)
                                if(k > 2*N_T)
                                    deltaV_min_3 = deltaV;
                                    k3_opt = k;
                                else
                                    deltaV_min_2 = deltaV;
                                    k2_opt = k;
                                end
                            else
                                deltaV_min_1 = deltaV;
                                k1_opt = k;
                            end
              
                            % Costo ottimale
                            deltaV_min = deltaV;

                            % Parametri orbitali del trasferimento
                            a_t_opt = a_t;
                            e_t_opt = e_t;
                            i_t_opt = i_t;
                            OM_t_opt = OM_t;
                            om_t_opt = om_t;
                            theta1_t_opt = theta1_t;
                            theta2_t_opt = theta2_t;

                            % Parametri orbitali di inizio e fine
                            % trasferimento
                            thetai_T = theta_T;
                            thetaf_A = theta_A_f;

                            % Tempo per il trasferimento
                            deltaT_t_opt = delta_T_t;

                            % Parametri per stabilire l'istante temporale
                            % del trasferimento ottimale
                            k_opt = k;
                            s_opt = s;
                            riv_T_opt = riv_T;

                        end
                    end
                end
            end
        end
    end
end

Time_3 = toc;
disp(['Tempo impiegato: ', num2str(Time_3), ' secondi']);

%% Parte grafica

delta_th = 0.005;
% Plot delle orbite eliocentriche
plotorbit(a_T, e_T, i_T, OM_T, om_T, th0_T, thf_T, delta_th, mu_s)
plotorbit(a_A, e_A, i_A, OM_A, om_A, th0_T, thf_T, delta_th, mu_s)

% Plot del Sole
scatter3(0, 0, 0, 100, 'b', 'filled');
hold on

% Plot del trasferimento interplanetario
plotorbit1(a_t_opt, e_t_opt, i_t_opt, OM_t_opt, om_t_opt, ...
    theta1_t_opt, theta2_t_opt, delta_th, mu_s)


[ri, ~ ] = par2car(a_t_opt, e_t_opt, i_t_opt, OM_t_opt, om_t_opt, theta1_t_opt, mu_s);
[rf, ~ ] = par2car(a_t_opt, e_t_opt, i_t_opt, OM_t_opt, om_t_opt, theta2_t_opt, mu_s);

scatter3(ri(1), ri(2), ri(3), 25, 'g', 'filled');
hold on
scatter3(rf(1), rf(2), rf(3), 25, 'r', 'filled');
hold on

legend("Orbita terrestre","Orbita dell'asteroide","Sole","Orbita di trasferimento ottimale", ...
    "Punto di partenza","Punto di arrivo")

%% Verifiche

mu_s = 1.327e11;

% a_t_opt = 1.848219275542337e+08 km
% e_t_opt = 0.185344860301687
% i_t_opt = 0.028484859392634 rad
% om_t_opt = 2.855993321445267 rad
% OM_t_opt = 1.231800840855091 rad
% theta1_t_opt = 0.282372259975278 rad
% theta2_t_opt = 4.752483335252830 rad
% deltaV_min_3 = 5.720441531561098 km/s (terza rivoluzione)
% deltaV_min_2 = 7.903437779646343 km/s (seconda rivoluzione)
% deltaV_min_1 = 7.989321300896556 km/s (prima rivoluzione)
% thetaf_A = 4.383695398283848 rad
% thetai_T = 2.602413028495237 rad
% k_opt = k_opt_3 = 7602
% k_opt_2 = 4799
% k_opt_1 = 1009

% a_t_opt = 1.848219275542337e+08;
% e_t_opt = 0.185344860301687;
% i_t_opt = 0.028484859392634;
% om_t_opt = 2.855993321445267;
% OM_t_opt = 1.231800840855091;
% theta1_t_opt = 0.282372259975278;
% theta2_t_opt = 4.752483335252830;
% 
% thetaf_A = 4.383695398283848;
% thetai_T = 2.602413028495237;

deltaT_t = TOF(a_t_opt, e_t_opt, theta1_t_opt, theta2_t_opt, mu_s);
% deltaT_t = 3.398569047730611e+07 s

% Determinazione delle velocità nei punti della manovra
[~ , vt_vect_1] = par2car(a_t_opt, e_t_opt, i_t_opt, OM_t_opt, om_t_opt, theta1_t_opt, mu_s);
[~ , vt_vect_2] = par2car(a_t_opt, e_t_opt, i_t_opt, OM_t_opt, om_t_opt, theta2_t_opt, mu_s);

[~ , v_vect_i] = par2car(a_T, e_T, i_T, OM_T, om_T, thetai_T, mu_s);
[~ , v_vect_f] = par2car(a_A, e_A, i_A, OM_A, om_A, thetaf_A, mu_s);

% Calcolo dei costi della manovra
deltaV1 = norm(vt_vect_1 - v_vect_i);
deltaV2 = norm(v_vect_f - vt_vect_2);
deltaV = deltaV1 + deltaV2;

% Dalla tolleranza temporale si può calcolare il massimo errore spaziale
% valutando l'angolo percorso nell'istante temporale e convertendolo nella
% lunghezza dell'arco di circonferenza associato

% % ------------------------------------------------------------------------
% % 1. Tolleranza temporale ammessa di 6 secondi (3 + 3)
% 
% thetaf_A = 4.326337702283748;
% 
% a_A = 1.966476025867008e+08;
% e_A = 0.283848775389335;
% mu_s = 1.327e11;
% 
% % Si parametrizza l'ellisse
% % x = a*cos(theta)
% % y = b*sin(theta)
% % Dunque l'arco sarà l'integrale rispetto a theta di sqrt(x^2 + y^2)
% 
% % Calcolo del semiasse minore dell'orbita dell'asteroide
% b_A = a_A * sqrt(1 - e_A^2);
% 
% thetaf_A_1 = thetaf_A;
% [~, thetaf_A_2] = KepleroDiretto(a_A, e_A, thetaf_A_1, 3, mu_s);
% 
% % Conversione da anomalia vera a anomalia eccentrica (E)
% E1 = 2 * atan2( sqrt(1 - e_A) * sin(thetaf_A_1 / 2), sqrt(1 + e_A) * cos(thetaf_A_1 / 2));
% E2 = 2 * atan2( sqrt(1 - e_A) * sin(thetaf_A_2 / 2), sqrt(1 + e_A) * cos(thetaf_A_2 / 2));
% 
% % Correzione per dominio dell'angolo
% if (E1 < 0)
%     E1 = E1 + 2*pi;
% end
% if (E2 < 0)
%     E2 = E2 + 2*pi;
% end
% 
% % Funzione integranda
% integranda = @(E) sqrt(a_A^2*sin(E).^2 + b_A^2*cos(E).^2);
% 
% % Calcolo dell'integrale numerico tra E1 ed E2
% s_1 = integral(integranda, E1, E2);  % Lunghezza in km
% 
% % Delta_s = 25.9797 km
% 
% % ------------------------------------------------------------------------
% 
% % 2. Tolleranza temporale ammessa di 4 secondi (2 + 2)
% 
% thetaf_A = 4.383695398283848;
% 
% % Si parametrizza l'ellisse
% % x = a*cos(theta)
% % y = b*sin(theta)
% % (Con theta che è l'anomalia eccentrica)
% % Dunque l'arco sarà l'integrale rispetto a theta di sqrt(x^2 + y^2)
% 
% % Calcolo del semiasse minore dell'orbita dell'asteroide
% b_A = a_A * sqrt(1 - e_A^2);
% 
% thetaf_A_1 = thetaf_A;
% [~, thetaf_A_2] = KepleroDiretto(a_A, e_A, thetaf_A_1, 2, mu_s);
% 
% % Conversione da anomalia vera a anomalia eccentrica (E)
% E1 = 2 * atan2( sqrt(1 - e_A) * sin(thetaf_A_1 / 2), sqrt(1 + e_A) * cos(thetaf_A_1 / 2));
% E2 = 2 * atan2( sqrt(1 - e_A) * sin(thetaf_A_2 / 2), sqrt(1 + e_A) * cos(thetaf_A_2 / 2));
% 
% % Correzione per dominio dell'angolo
% if (E1 < 0)
%     E1 = E1 + 2*pi;
% end
% if (E2 < 0)
%     E2 = E2 + 2*pi;
% end
% 
% % Funzione integranda
% integranda = @(E) sqrt(a_A^2*sin(E).^2 + b_A^2*cos(E).^2);
% 
% % Calcolo dell'integrale numerico tra E1 ed E2
% s_2 = abs(integral(integranda, E1, E2));  % Lunghezza in km
% 
% % Delta_s = 0.31597 km

%% ------------------------------------------------------------------------

% 3. Tolleranza temporale ammessa di 2 secondi (1 + 1)

thetaf_A = 4.383695398283848;

% Si parametrizza l'ellisse
% x = a*cos(theta)
% y = b*sin(theta)
% Dunque l'arco sarà l'integrale rispetto a theta di sqrt(x^2 + y^2)

% Calcolo del semiasse minore dell'orbita dell'asteroide
b_A = a_A * sqrt(1 - e_A^2);

thetaf_A_1 = thetaf_A;
[~, thetaf_A_2] = KepleroDiretto(a_A, e_A, thetaf_A_1, 1, mu_s);

% Conversione da anomalia vera a anomalia eccentrica (E)
E1 = 2 * atan2( sqrt(1 - e_A) * sin(thetaf_A_1 / 2), sqrt(1 + e_A) * cos(thetaf_A_1 / 2));
E2 = 2 * atan2( sqrt(1 - e_A) * sin(thetaf_A_2 / 2), sqrt(1 + e_A) * cos(thetaf_A_2 / 2));

% Correzione per dominio dell'angolo
if (E1 < 0)
    E1 = E1 + 2*pi;
end
if (E2 < 0)
    E2 = E2 + 2*pi;
end

% Funzione integranda
integranda = @(E) sqrt(a_A^2*sin(E).^2 + b_A^2*cos(E).^2);

% Calcolo dell'integrale numerico tra E1 ed E2
s_3 = abs(integral(integranda, E1, E2));  % Lunghezza in km

% Instabilità della soluzione
% ------------------------------------------------------------------------


% ------------------------------------------------------------------------
% Rappresentazione grafica dell'ottimizzazione dei costi

% Mashgrid dei parametri
[TH1, W] = meshgrid(theta_T_vect, om_t_vect);

% ------------------------------------------------------------------------

% Prima rivoluzione terrestre
figure
surf(TH1, W, deltaV_map_1',  'EdgeColor', 'none')
xlabel('Anomalia vera terrestre iniziale [rad]')
ylabel('Anomalia di perigeo orbita di trasferimento [rad]')
zlabel('\DeltaV [m/s]')
title('Andamento del costo \DeltaV in funzione di \theta_i e di \omega_t')
colorbar;
view(135, 30) % Per vedere da una prospettiva gradevole

figure;
contourf(TH1, W, deltaV_map_1', 50);  % 50 linee di livello
xlabel('Anomalia vera terrestre iniziale [rad]');
ylabel('Anomalia di perigeo orbita di trasferimento [rad]');
title('Contorni del costo \DeltaV');
colorbar;
hold on;

% Evidenzia punto ottimo
scatter(thetai_T, om_t_opt, 100, 'r', 'filled');
text(thetai_T, om_t_opt, ' \leftarrow Ottimo', 'Color', 'red');

% ------------------------------------------------------------------------

% Seconda Rivoluzione terrestre
figure
surf(TH1, W, deltaV_map_2',  'EdgeColor', 'none')
xlabel('Anomalia vera terrestre iniziale [rad]')
ylabel('Anomalia di perigeo orbita di trasferimento [rad]')
zlabel('\DeltaV [m/s]')
title('Andamento del costo \DeltaV in funzione di \theta_i e di \omega_t')
colorbar;
view(135, 30) % Per vedere da una prospettiva gradevole

figure;
contourf(TH1, W, deltaV_map_2', 50);  % 50 linee di livello
xlabel('Anomalia vera terrestre iniziale [rad]');
ylabel('Anomalia di perigeo orbita di trasferimento [rad]');
title('Contorni del costo \DeltaV');
colorbar;
hold on;

% Evidenzia punto ottimo
scatter(thetai_T, om_t_opt, 100, 'r', 'filled');
text(thetai_T, om_t_opt, ' \leftarrow Ottimo', 'Color', 'red');

% ------------------------------------------------------------------------

% Terza Rivoluzione terrestre
figure
surf(TH1, W, deltaV_map_3',  'EdgeColor', 'none')
xlabel('Anomalia vera terrestre iniziale [rad]')
ylabel('Anomalia di perigeo orbita di trasferimento [rad]')
zlabel('\DeltaV [m/s]')
title('Andamento del costo \DeltaV in funzione di \theta_i e di \omega_t')
colorbar;
view(135, 30) % Per vedere da una prospettiva gradevole

figure;
contourf(TH1, W, deltaV_map_3', 50);  % 50 linee di livello
xlabel('Anomalia vera terrestre iniziale [rad]');
ylabel('Anomalia di perigeo orbita di trasferimento [rad]');
title('Contorni del costo \DeltaV');
colorbar;
hold on;

% Evidenzia punto ottimo
scatter(thetai_T, om_t_opt, 100, 'r', 'filled');
text(thetai_T, om_t_opt, ' \leftarrow Ottimo', 'Color', 'red');
% ------------------------------------------------------------------------

figure('Position', [100, 100, 1200, 800]);
    
% 3D orbit plot
subplot(2, 2, [1, 2]);
plotorbit(a_T, e_T, i_T, OM_T, om_T, th0_T, thf_T, delta_th, mu_s);
hold on;
plotorbit(a_A, e_A, i_A, OM_A, om_A, th0_T, thf_T, delta_th, mu_s);
plotorbit1(a_t_opt, e_t_opt, i_t_opt, OM_t_opt, om_t_opt, ...
    theta1_t_opt, theta2_t_opt, delta_th, mu_s);

% Plot Sun and transfer points
scatter3(0, 0, 0, 200, 'yellow', 'filled', 'MarkerEdgeColor', 'black');
[ri, ~] = par2car(a_t_opt, e_t_opt, i_t_opt, OM_t_opt, om_t_opt, theta1_t_opt, mu_s);
[rf, ~] = par2car(a_t_opt, e_t_opt, i_t_opt, OM_t_opt, om_t_opt, theta2_t_opt, mu_s);
scatter3(ri(1), ri(2), ri(3), 100, 'green', 'filled');
scatter3(rf(1), rf(2), rf(3), 100, 'red', 'filled');

legend('Earth Orbit', 'Asteroid Orbit', 'Transfer Orbit', 'Sun', 'Departure', 'Arrival', ...
    'Location', 'best');
title('Optimal Interplanetary Transfer (Parallelized)');
grid on;
axis equal;

% Cost surface plot
subplot(2, 2, 3);
[TH1, W] = meshgrid(theta_T_vect, om_t_vect);
deltaV_map_plot = deltaV_map_1';
deltaV_map_plot(deltaV_map_plot == inf) = NaN; % Replace inf with NaN for plotting

surf(TH1, W, deltaV_map_plot, 'EdgeColor', 'none');
xlabel('Earth True Anomaly [rad]');
ylabel('Transfer Periapsis Argument [rad]');
zlabel('\DeltaV [m/s]');
title('Cost Surface');
colorbar;
view(45, 30);

% Contour plot
subplot(2, 2, 4);
contourf(TH1, W, deltaV_map_plot, 20);
xlabel('Earth True Anomaly [rad]');
ylabel('Transfer Periapsis Argument [rad]');
title('Cost Contours');
colorbar;
hold on;

% Mark optimum
scatter(thetai_T, om_t_opt, 100, 'red', 'filled', 'MarkerEdgeColor', 'white');
text(thetai_T, om_t_opt, sprintf('  Optimum\n  ΔV=%.0f m/s', deltaV_min), ...
    'Color', 'red', 'FontWeight', 'bold');
