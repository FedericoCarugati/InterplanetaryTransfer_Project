clc
clear
close all
% figure

% Il secondo scenario si propone di ottimizzare un trasferimento diretto da
% un'orbita eliocentrica della Terra ad un'orbita eliocentrica di un
% asteroide assegnato 

% Ottimizzazione dei tempi

% Tempo stimato per lo script: 1-2 minuti

% Tempo totale del trasferimento : 0.6181 giorni
% Costo propulsivo (DeltaV) : 106.69 km/s

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

th0 = 0;
thf = 2*pi;
delta_th = 0.005;

%% Parametri orbitali dell'orbita eliocentrica dell'asteroide

% Semiasse maggiore
a_A =  1.324170*UA;
% Eccentricità
e_A = 0.280134;
% Inclinazione
i_A =  1.62*(pi/180);
% RAAN - Ascensione retta del nodo ascendente
OM_A = 69.08*(pi/180);
% Anomalia del pericentro
om_A = 162.84*(pi/180);

th0 = 0;
thf = 2*pi;
delta_th = 0.005;

%% Parametri orbitali dell'orbita di trasferimento interplanetario

n = 100;

% Il piano orbitale dell'orbita di trasferimento è dipendente dal punto di
% partenza sull'orbita terrestre e quello di arrivo sull'orbita
% dell'asteroide (theta1_i e theta2_f)

% Anomalia vera iniziale sull'orbita eliocentrica terrestre
theta1_i_vect = linspace(0, 2*pi, n);
% Anomalia vera finale sull'orbita eliocentrica dell'asteroide
theta2_f_vect = linspace(0, 2*pi, n);

% Anomalia vera dell'orbita di trasferimento
om_t_vect = linspace(0, 2*pi, n); 

% Si procede con un'ottimizzazione dei costi con dei cicli for annidati che
% facciano variare i parametri delle anomalie vere (iniziale e finale) e
% dell'anomalia di pericentro dell'orbita di trasferimento

% Inizializzazione della variabile per l'ottimizzazione dei costi
deltaT_min = 0;
deltaT_vect = [];

tic;

for(j = 1:n)
    theta1_i = theta1_i_vect(j);
    for(k = 1:n)
        theta2_f = theta2_f_vect(k);
        [r_vect_i, v_vect_i] = par2car(a_T, e_T, i_T, OM_T, om_T, theta1_i, mu_s);
        [r_vect_f, v_vect_f] = par2car(a_A, e_A, i_A, OM_A, om_A, theta2_f, mu_s);

        % Si definisce il piano orbitale di trasferimento che congiunge i
        % due punti dello spazio definiti dalle posizioni iniziale e finale
        % sulle due orbite (date dalle anomalie)
        h_t_vect = cross(r_vect_i, r_vect_f)/norm(cross(r_vect_i, r_vect_f));

        i_vect = [1;0;0];
        j_vect = [0;1;0];
        k_vect = [0;0;1];

        % Inclinazione dell'orbita
        i_t = acos(dot(h_t_vect,k_vect));

        % Determinazione della linea dei nodi
        N_vect = cross(k_vect,h_t_vect)/norm(cross(k_vect,h_t_vect));

        % Calcolo dell'ascensione retta del nodo ascenendente 
        if(dot(N_vect,j_vect)>=0)
            OM_t = acos(dot(N_vect,i_vect));
        else
            OM_t = 2*pi - acos(dot(N_vect,i_vect));
        end

        for(s = 1:n)
            om_t = om_t_vect(s);

            % Si possono calcolare tutti i parametri dell'orbita di
            % trasferimento considerando i costi della missione

            % Matrici di rotazione

            R_om = [cos(OM_t)    sin(OM_t)     0;
                    -sin(OM_t)   cos(OM_t)     0;
                        0           0          1];

            R_i = [     1           0          0;
                        0        cos(i_t)   sin(i_t);
                        0       -sin(i_t)   cos(i_t)];

            R_w = [cos(om_t_vect(s)) sin(om_t_vect(s))           0;
                  -sin(om_t_vect(s)) cos(om_t_vect(s))           0;
                          0                  0                   1];

            T = R_w*R_i*R_om;

            % Definizione dei vettori posizione dei punti 1 e 2 sull'orbita 
            % di trasferimnto nel sistema perifocale
            rt_1_vect = T*r_vect_i;
            rt_2_vect = T*r_vect_f;

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

            if((e_t > 0) && (1 > e_t))
                % Semiasse maggiore
                a_t = norm(rt_1_vect)*((1 + e_t*cos(theta1_t))/(1 - e_t^2));
                % Raggio di pericentro 
                rp_t = a_t*(1 - e_t);
           
                if(rp_t > 0.7*UA)
                    % Determinazione delle velocità nei punti della manovra
                    [~ , vt_vect_1] = par2car(a_t, e_t, i_t, OM_t, om_t, theta1_t, mu_s);
                    [~ , vt_vect_2] = par2car(a_t, e_t, i_t, OM_t, om_t, theta2_t, mu_s);
        
                    % Durata del trasferimento
                    deltaT = TOF(a_t, e_t, theta1_t, theta2_t, mu_s);

                    % Calcolo dei costi della manovra
                    deltaV1 = norm(vt_vect_1 - v_vect_i);
                    deltaV2 = norm(v_vect_f - vt_vect_2);
                    deltaV = deltaV1 + deltaV2;
                    deltaT_vect = [deltaT_vect deltaV];
        
                    % Selezione del trasferimento ottimo
                    if(deltaT_min == 0 || deltaT_min > deltaT)
                        deltaT_min = deltaT;
                        deltaV_Tmin = deltaV;
                        a_t_opt = a_t;
                        e_t_opt = e_t;
                        i_t_opt = i_t;
                        OM_t_opt = OM_t;
                        om_t_opt = om_t;
                        theta1_t_opt = theta1_t;
                        theta2_t_opt = theta2_t;
                        k_opt = k;
                        j_opt = j;
                        s_opt = s;
                        theta1_i_opt = theta1_i_vect(j);
                    end
                end
            end
        end
    end
end

Time = toc;
disp(['Tempo impiegato: ', num2str(Time), ' secondi']);

S = [0.7500 0.3250 0.0980];

% Plot delle orbite
plotorbit(a_T, e_T, i_T, OM_T, om_T, th0, thf, delta_th, mu_s)
plotorbit(a_A, e_A, i_A, OM_A, om_A, th0, thf, delta_th, mu_s)

% Plot del Sole
scatter3(0, 0, 0, 100, S,'filled');
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

% Da tenere presente che i valori delle anomalie vere di inizio e fine 
% trasferimento possono condizionare notevolmente i costi degli scenari
% precedenti, andando ad intaccare il costo complessivo della missione

legend("Earth Orbit", "Asteroid Orbit", "Sun", "Optimal Transfer Orbit", ...
       "Departure Point", "Arrival Point")

% deltaT_min = 0.6181 giorni
% deltaV_Tmin = 106.69 km/s


