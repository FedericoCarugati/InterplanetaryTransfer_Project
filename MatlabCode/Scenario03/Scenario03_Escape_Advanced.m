clc
close all
clear
% figure

% Calcolo della manovra di fuga ottimale
% deltaV_fuga = 2.936408513453976 km/s

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

% Raggio della Terra
R_T = 6373;

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

p_f = (1 - e_f^2)*a_f; % Semilato retto

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


%% INIZIO DELLO SCENARIO 3

%% Calcolo dei parametri dell'Orbita Iperbolica di fuga

% In questo caso non si considera più la semplificazione di orbita di
% parcheggio complanare all'orbita iperbolica di fuga dunque viene anche a
% decadere l'ipotesi che il perigeo delle due orbite coincida

% In questo caso la necessità è quella di far partire l'orbita di fuga in
% modo che raggiunga la v_inf con la direzione coerente a quanto calcolato
% nello scenario 2

% Calcolo del vettore v_inf_vect nel sistema di riferimento geocentrico
% partendo dal calcolo della velocità relativa nel sistema eliocentrico
% eclittico e poi si applica una rotazione epsilon per il cambio di piano
% (cambio del sistema di riferimento)

v_inf_vect_eclipt = vt_vect_1 - vT_vect_i;

epsilon = 23.45*(pi/180); % Angolo di inclinazione del sistema terrestre ù

% Matrice di rotazione:
T = [1          0             0;
     0    cos(epsilon)   sin(epsilon);
     0   -sin(epsilon)   cos(epsilon)];

% Velocità infinito in ECI
v_inf_vect = T'*v_inf_vect_eclipt; 

% Dunque, si calcola la direzione di uscita della velocità
r_inf_vers = v_inf_vect/norm(v_inf_vect);

% Questa direzione sarà la direzione dell'asindoto riferito all'iperbole di
% fuga (sarà una retta definita nello spazio):
% r = @(x,y) i_vect_inf(3)/i_vect_inf(2)*y + i_vect_inf(3)/i_vect_inf(1)*x;

% Si deve ricavare un'iperbole che interseca l'ellisse e che riesca a
% soddisfare l'ipotesi della direzione della velocità all'infinito

v_inf = norm(v_inf_vect); % Norma della v_inf
a_h = - mu/v_inf^2; % Semiasse maggiore dell'iperbole

% A tale scopo si fa ciclare un vettore di posizioni da cui far partire
% l'iperbole e si valuta l'esistenza di un'iperbole che, partendo da quel
% punto, soddisfi le richieste iniziali e che sia anche ottimale a livelli
% di costo

% Definiamo l'array di vettori che identificano le intersezioni tra
% l'iperbole e l'orbita di parcheggio
N = 500;
r_h_vect = zeros(3,N); % Inizializzazione del vettore posizione di intersezione
v_e_vect = zeros(3,N);
r_h_vers = zeros(3,N); % Inizializzazione del versore per la direzione di intersezione
theta_f_vect = linspace(0, 2*pi, N); % Inizializzazione del vettore delle anomalie vere

% Calcolo dell'array di vettori posizione
for(i=1:N)
    [r_h_vect(:,i), v_e_vect(:,i)] = par2car(a_f, e_f, i_f, OM_f, om_f, theta_f_vect(i), mu);
    r_h_vers(:,i) = r_h_vect(:,i)/norm(r_h_vect(:,i));
end

% it_vect = [];
% x_vect_v = [];
x_sol_par_vect = [];
deltaV_min = 0;
t = 0;
deltaV_vect = [];
r_h_v_poss = [];

% Ciclo di iterazioni per l'ottimizzazione

options = optimoptions('fsolve', 'Display', 'off');

for(i=1:N) 
    r_h = norm(r_h_vect(:,i)); % Modulo del vettore posizione di intersezione
    alpha = acos(dot(r_h_vers(:,i), r_inf_vers)); 
    if(0 > dot(r_h_vers(:,i), r_inf_vers))
        alpha = 2*pi - alpha;
    end

    % ------------------------------------------------------------------
    % METODO CON LA FUNZIONE DI NEWTON
    % Definizione del vettore funzione
    % F1 = @(e, theta_h, theta_inf) [r_h*(1 + e*cos(theta_h)) - a_h*(1 - e^2);
    %                                cos(theta_inf) + 1/e;
    %                                theta_inf - alpha - theta_h];
    % % Jacobiano
    % J1 = @(e, theta_h, theta_inf) [(r_h*cos(theta_h) + 2*a_h*e)  (-r_h*e*sin(theta_h))          0;
    %                                           -1/e^2                      0              (-sin(theta_inf));
    %                                             0                        -1                     1];
    % 
    % % Punto di partenza
    % x0 = 1.3*ones(3);
    % % Calcolo dei parametri
    % [x_vect, it] = newton_sys3(x0, 1e6, 1e-8, F1, J1);
    % 
    % % Aggiornamenti del vettore delle soluzioni
    % it_vect = [it_vect, it];
    % x_vect_v = [x_vect_v x_vect(:,it)];
    % ------------------------------------------------------------------

    % METODO CON LA FUNZION fsolve DI MATLAB
    % Definizione del vettore funzione
    F1 = @(x) [r_h*(1 + x(1)*cos(x(2))) - a_h*(1 - x(1)^2);
               cos(x(3)) + 1/x(1);
               x(3) - alpha - x(2)];
    x0 = ones(3,1);

    x_sol_par = fsolve(F1, x0, options);
    x_sol_par_vect = [x_sol_par_vect x_sol_par];

    % Parametri dell'iperbole
    e_h = x_sol_par(1); % Eccentricità dell'iperbole
    theta_h_i = x_sol_par(2); % Anomalia sull'iperbole riferita al punto di intersezione
    theta_inf_i = x_sol_par(3); % Anomalia vera all'infinito dell'iperbole    

    % Data l'eccentricità dell'iperbole si calcola il raggio di pericentro
    rp_h = a_h*(1 - e_h);

    if(theta_inf_i > pi/2 && pi > theta_inf_i && rp_h > R_T)
        % Si cerca di ricavare i restanti parametri dell'orbita iperbolica

        % Si calcola il modulo del vettore del momento angolare specifico e
        % si può anche ricavare la direzione tramite una normalizzazione
        % del prodotto vettoriale tra il vettore velocità v_inf_vect e il
        % vettore posizione nel punto di intersezione tra le orbite
        p_h = a_h*(1 - e_h^2); % Semilato retto
        h_h = sqrt(mu*p_h); % Modulo del momento angolare specifico
        h_h_vers = cross(r_h_vect(:,i), v_inf_vect)/norm(cross(r_h_vect(:,i), v_inf_vect));
        h_h_vect = h_h*h_h_vers; % Vettore momento angolare

        % Calcolo alternativo dell'inclinazione del piano orbitale
        K = [0;0;1];
        i_h = acos(dot(h_h_vect, K)/norm(h_h_vect));
        % Questo angolo riferito all'inclinazione dell'orbita dovrebbe
        % essere approssimativamente uguale a quello calcolato
        % precedentemente

        % Nodo ascendente dell'iperbole
        N_h = cross(K, h_h_vect)/norm(cross(K, h_h_vect));
        I = [1;0;0];
        % Ascensione retta del nodo ascendente
        OM_h = acos(dot(N_h, I));
        j = [0;1;0];
        if(0 > dot(N_h, j))
            OM_h = 2*pi - OM_h;
        end

        % Calcolo dell'anomalia di pericentro dell'iperbole:

        % PRIMO METODO:
        beta = asin(dot(cross(N_h, r_h_vers(:,i)), h_h_vers));
        if(0 > dot(N_h, r_h_vers(:,i)))
            beta = pi - beta;
        end
        om_h = beta - theta_h_i; % Anomalia di pericentro

        % SECONDO METODO:

        % ------------------------------------------------------------------
        % Metodo implementato tramite la funzione di Newton per la
        % risoluzione del sistema non lineare
        % F2 = @(e1, e2, e3) [r_h_vers(1,i)*e1 + r_h_vers(2,i)*e2 + r_h_vers(3,i)*e3 - cos(theta_h_i);
        %                     r_inf_vers(1)*e1 + r_inf_vers(2)*e2 + r_inf_vers(3)*e3 - cos(theta_inf_i);
        %                     e1^2 + e2^2 + e3^2 - 1;];
        % % Jacobiano
        % J2 = @(e1, e2, e3) [r_h_vers(1,i)   r_h_vers(2,i)   r_h_vers(3,i);
        %                     r_inf_vers(1)   r_inf_vers(2)   r_inf_vers(3);
        %                        2*e1             2*e2            2*e3];
        % 
        % % Punto di partenza
        % x0 = 1.5*ones(3);
        % % Calcolo dei coefficienti del vettore eccentricità
        % [x_vect_2, it_2] = newton_sys3(x0, 1e3, 1e-5, F2, J2);
        % % Costruzione del versore eccentricità
        % e_h_vers = [x_vect_2(1,it_2); x_vect_2(2,it_2); x_vect_2(3,it_2)];
        % ------------------------------------------------------------------

        % Metodo implementato tramite la funzione fsolve di MATLAB per la
        % risoluzione del sistema non lineare
        F2 = @(x) [r_h_vers(1,i)*x(1) + r_h_vers(2,i)*x(2) + r_h_vers(3,i)*x(3) - cos(theta_h_i);
                   r_inf_vers(1)*x(1) + r_inf_vers(2)*x(2) + r_inf_vers(3)*x(3) - cos(theta_inf_i);
                   x(1)^2 + x(2)^2 + x(3)^2 - 1;];

        x_sol_e = fsolve(F2, x0, options);

        % Costruzione del versore eccentricità
        e_h_vers = [x_sol_e(1); x_sol_e(2); x_sol_e(3)];

        % Questa anomalia vera deve coincidere con quella già calcolata ai
        % passi precedenti
        theta_h_1 = acos(dot(e_h_vers, r_h_vers(:,i)));
        if(0 > dot(r_h_vect(:,i), v_e_vect(:,i)))
            theta_h_1 = 2*pi - theta_h_1;
        end

        % Si verifica che i due metodi per ricavare l'anomalia dal perigeo siano
        % uguali
        om_h_2 = acos(dot(N_h, e_h_vers));
        if(0 > dot(e_h_vers, K))
            om_h_2 = 2*pi - om_h_2;
        end
        % Si verifica che om_h_2 = om_h (i due metodi sono identici)

        % Si calcola il vettore velocità dell'orbita di fuga nel punto di
        % intersezione tra le orbite

        % Matrici di rotazione:

        R_OM = [ cos(OM_h)    sin(OM_h)       0;
                -sin(OM_h)    cos(OM_h)       0;
                    0            0            1];  

        R_i = [     1            0              0;
                    0         cos(i_h)      sin(i_h);
                    0        -sin(i_h)      cos(i_h)]; 

        R_om = [ cos(om_h)     sin(om_h)        0;
                -sin(om_h)     cos(om_h)        0;
                    0             0             1];    

        T = R_OM*R_i*R_om;     

        % Vettore velocità nel sistema perifocale
        v_pf = sqrt(mu/p_h)*[-sin(theta_h_1); e_h + cos(theta_h_1); 0];        
        % Vettori velocità nel sistema di rif. celeste
        v_h_vect = T'*v_pf;

        % Dunque si calcola il costo della manovra
        deltaV = norm(v_h_vect - v_e_vect(:,i));
        deltaV_vect = [deltaV_vect deltaV];

        r_h_v_poss = [r_h_v_poss r_h];

        if(deltaV_min > deltaV || t==0)
            t = t+1;

            i_opt = i; % Iterazione ottimale

            % Aggiornamento dei parametri ottimali

            deltaV_min = deltaV; % Costo ottimale

            % Parametri iperbole di fuga
            e_h_opt = e_h; % Eccentricità
            e_h_vers_opt = e_h_vers; % Direzione del vettore eccentricità
            om_h_opt = om_h; % Anomalia di pericentro
            OM_h_opt = OM_h; % Ascensione retta del nodo ascendente
            i_h_opt = i_h; % Inclinazione del piano orbitale

            r_h_opt = r_h;

            % Anomalie delle manovre impulsive
            theta_h_opt = theta_h_i; % Anomalia vera di inizio fuga sull'orbita iperbolica
            theta_h_1_opt = theta_h_1; % Uguale a theta_h_opt
            theta_inf_opt = theta_inf_i; % Anomalia vera all'infinito
            alpha_opt = alpha; % Angolo tra il vettore posizione iniziale e il versore di uscita
            theta_f_i_opt = theta_f_vect(i); % Anomali vera di inizio fuga sull'orbita ellittica
        end
    end
end

deltaV_fuga = deltaV_min; % Costo della fuga

% Verifica che sia andato tutto per il verso giusto
r_e_f_vect = par2car(a_f, e_f, i_f, OM_f, om_f, theta_f_i_opt, mu);
r_h_i_vect = par2car(a_h, e_h_opt, i_h_opt, OM_h_opt, om_h_opt, theta_h_opt, mu);

deltaV_TOT = deltaV_tot + deltaV_fuga

%% Parte grafica
% Disegno del pianeta
scatter3(0, 0, 0, 80, 'blue', 'filled');
hold on

% Disegno del punto
scatter3(r_i(1), r_i(2), r_i(3), 25, 'g', 'filled');
hold on

p_h = a_h*(1 - e_h_opt^2); 
% thf_h_1 = theta_inf_opt;
thf_h_1 = acos((p_h/(3*a_f) - 1)/e_h_opt);
th0_h_1 = - thf_h_1;

% Disegno dell'orbita iniziale e finale
plotorbit(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi, delta_th, mu)
plotorbit(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi, delta_th, mu)

% Disegno dell'orbita di fuga
plotorbit_hyperbolic(a_h, e_h_opt, i_h_opt, OM_h_opt, om_h_opt, th0_h_1, thf_h_1, delta_th, mu)

% Animazione dei trasferimenti
plotorbit_improved(a_i, e_i, i_i, OM_i, om_i, theta_i, pi, delta_th, mu)
plotorbit_improved(ra_i, 0, i_i, OM_i, 0, (om_i + pi), theta_circ_1_ell, delta_th, mu)
plotorbit_improved(a_t, e_t, i_i, OM_i, theta_circ_1_ell, 0, pi, delta_th, mu)
plotorbit_improved(ra_f, 0, i_f, OM_f, 0, theta_circ_2, (om_f + pi), delta_th, mu)
plotorbit_improved(a_f, e_f, i_f, OM_f, om_f, pi, theta_f_i_opt, delta_th, mu)

thf_h = acos((p_h/(2.5*a_f) - 1)/e_h_opt);

th0_i = 0;
thf_i = 2*pi;

% Animazione della manovra di fuga
plotorbit_improved_hyperbolic(a_h, e_h_opt, i_h_opt, OM_h_opt, om_h_opt, theta_h_1_opt, thf_h, delta_th, mu)

legend('Earth', 'Starting point', 'Initial orbit', 'Parking orbit', ...
    'Hyperbolic escape orbit', 'Transfer from initial point to start of circularization', ...
    'First impulsive maneuver', 'First circular orbit', 'Second impulsive maneuver', ...
    'Connecting ellipse between circular orbits', 'Third impulsive maneuver', ...
    'Second circular orbit', 'Fourth impulsive maneuver', 'Transfer on final orbit', ...
    'Escape maneuver', 'Hyperbolic escape trajectory')

% n = length(deltaV_vect);
% s = 0;
% for(i=1:(n-1))
%     if(r_h_v_poss(i) > r_h_v_poss(i+1) && s==0)
%         s = 1;
%         n_stop = i;
%     end
% end
% 
% deltaV_vect1 = deltaV_vect(1:n_stop);
% r_h_v_poss1 = r_h_v_poss(1:n_stop);
% deltaV_vect2 = deltaV_vect(n_stop+1:end);
% r_h_v_poss2 = r_h_v_poss(n_stop+1:end);

% figure 
% plot(r_h_v_poss1, deltaV_vect1)
% xlabel('Distanza del punto di inizio fuga da attrattore (r_h)')
% ylabel('Costo manovra (deltaV)')
% hold on
% plot(r_h_v_poss2, deltaV_vect2)
% scatter(r_h_opt, deltaV_vect(i_opt))
% legend('Prima metà orbita di parcheggio','Seconda metà orbita di parcheggio',...
%     'Ottimo')

