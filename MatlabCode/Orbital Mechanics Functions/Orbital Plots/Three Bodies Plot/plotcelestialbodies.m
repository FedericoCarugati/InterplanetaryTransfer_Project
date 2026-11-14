function plotcelestialbodies(a_1, e_1, i_1, OM_1, om_1, th0_1, thf_1, ...
    a_2, e_2, i_2, OM_2, om_2, th0_2, dt, mu, theta1_sat, theta2_sat, ...
    thetaf_2, thetai_1, a_s, e_s, i_s, om_s, OM_s)

% Visualizzazione animata dell'orbita in 3D con una scia persistente
%
% ------------------------------------------------------------------
% Argomenti di Input
% a_1              Semi-asse maggiore                       [km]
% e_1              Eccentricità                             [-]
% i_1              Inclinazione                             [rad]
% OM_1             Ascensione retta del nodo ascendente     [rad]
% om_1             Anomalia del pericentro                  [rad]
% th0_1            Anomalia vera iniziale                   [rad]
% thf_1            Anomalia vera finale                     [rad]
% a_2              Semi-asse maggiore                       [km]
% e_2              Eccentricità                             [-]
% i_2              Inclinazione                             [rad]
% OM_2             Ascensione retta del nodo ascendente     [rad]
% om_2             Anomalia del pericentro                  [rad]
% th0_2            Anomalia vera iniziale                   [rad]
% dt               Passo temporale                          [rad]
% mu               Costante gravitazionale                  [km^3/s^2]
% theta1_sat       Anomalia vera iniziale satellite         [rad]
% theta2_sat       Anomalia vera finale satellite           [rad]
% thetaf_2         Anomalia vera finale manovra             [rad]
% thetai_1         Anomalia vera iniziale manovra           [rad]
% a_s              Semi-asse maggiore                       [km]
% e_s              Eccentricità                             [-]
% i_s              Inclinazione                             [rad]
% OM_s             Ascensione retta del nodo ascendente     [rad]
% om_s             Anomalia del pericentro                  [rad]
%
% ------------------------------------------------------------------

tic;
% Controllo sugli input
if a_1 <= 0 || e_1 < 0 || e_1 >= 1 || mu <= 0
    error('Parametri orbitali non validi.');
end

% Controllo sugli input
if a_2 <= 0 || e_2 < 0 || e_2 >= 1 || mu <= 0
    error('Parametri orbitali non validi.');
end

% Tempo per arrivare all'anomalia iniziale dal perigeo
Ti_1 = TOF(a_1, e_1, 0, th0_1, mu);
% Ti_2 = TOF(a_2, e_2, 0, th0_2, mu);
% Tempo orbitale dall'anomalia iniziale 
T_1 = TOF(a_1, e_1, th0_1, thf_1, mu);
% T_2 = TOF(a_2, e_2, th0_2, th0_2 + 2*pi, mu);
% Tempo complessivo dal perigeo
Tf_1 = T_1 + Ti_1;
% Tf_2 = T_2 + Ti_2;

T2_tot = TOF(a_2, e_2, 0, 2*pi, mu);
N_tot = round(T2_tot/dt);

% Vettore dei tempi
t_vect = Ti_1 : dt : Tf_1; % Si considerano gli stessi tempi delle orbite
% t_vect = Ti_2 : dt : Tf_2;

T_t = TOF(a_s, e_s, theta1_sat, theta2_sat, mu);
t_vect_2 = 0 : dt : T_t;

N = length(t_vect);
n = length(t_vect_2);

% Inizializzazione del vettore delle anomalie vere
theta_1 = zeros(1, N);
theta_2 = zeros(1, N);
theta_s = zeros(1, n);
theta_1(1) = th0_1;
theta_2(1) = th0_2;
theta_s(1) = theta1_sat;

% Calcolo posizione ad istanti temporali regolari 

r_1 = zeros(3, N);
v_1 = zeros(3, N);
[r_i, v_i] = par2car(a_1, e_1, i_1, OM_1, om_1, th0_1, mu);
r_1(:, 1) = r_i;
v_1(:, 1) = v_i;

r_2 = zeros(3, N);
v_2 = zeros(3, N);
[r_i, v_i] = par2car(a_2, e_2, i_2, OM_2, om_2, th0_2, mu);
r_2(:, 1) = r_i;
v_2(:, 1) = v_i;

r_s = zeros(3, n);
v_s = zeros(3, n);
[r_i, v_i] = par2car(a_s, e_s, i_s, OM_s, om_s, theta1_sat, mu);
r_s(:, 1) = r_i;
v_s(:, 1) = v_i;

for j = 1:N-1
    % Calcolo dell'anomalia vera successiva
    [~, theta2] = KepleroDiretto(a_1, e_1, theta_1(j), dt, mu);
    theta_1(j+1) = theta2;
    % Calcolo del vettore corrispondente
    [r_i, v_i] = par2car(a_1, e_1, i_1, OM_1, om_1, theta_1(j+1), mu);
    r_1(:, j+1) = r_i;
    v_1(:, j+1) = v_i;

    % Calcolo dell'anomalia vera successiva
    [~, theta2] = KepleroDiretto(a_2, e_2, theta_2(j), dt, mu);
    theta_2(j+1) = theta2;
    % Calcolo del vettore corrispondente
    [r_i, v_i] = par2car(a_2, e_2, i_2, OM_2, om_2, theta_2(j+1), mu);
    r_2(:, j+1) = r_i;
    v_2(:, j+1) = v_i;
    
end

for j = 1:n
    % Calcolo dell'anomalia vera successiva
    [~, theta2] = KepleroDiretto(a_s, e_s, theta_s(j), dt, mu);
    theta_s(j+1) = theta2;
    % Calcolo del vettore corrispondente
    [r_i, v_i] = par2car(a_s, e_s, i_s, OM_s, om_s, theta_s(j+1), mu);
    r_s(:, j+1) = r_i;
    v_s(:, j+1) = v_i;
end

elapsedTime = toc;
disp(['Tempo impiegato: ', num2str(elapsedTime), ' secondi']);

% Creazione della figura
hold on
grid on
axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('3D Animation')
view(3)
rotate3d on

% Punti inziali Terra e asteroide
scatter3(r_1(1,1), r_1(2,1), r_1(3,1), 25, '*', 'g')
hold on
scatter3(r_2(1,1), r_2(2,1), r_2(3,1), 25, '*', 'r')
hold on

% Inizializzazione delle traiettorie
traj_1 = plot3(r_1(1,1), r_1(2,1), r_1(3,1), 'LineWidth', 0.5);
traj_2 = plot3(r_2(1,1), r_2(2,1), r_2(3,1), 'LineWidth', 0.5);
sat = plot3(r_s(1,1), r_s(2,1), r_s(3,1), 'LineWidth', 0.5);

% Corpi celesti
corpo_1 = scatter3(r_1(1,1), r_1(2,1), r_1(3,1), 40, 'g', 'filled');
corpo_2 = scatter3(r_2(1,1), r_2(2,1), r_2(3,1), 40, 'r', 'filled');
satellite = scatter3(r_s(1,1), r_s(2,1), r_s(3,1), 40, 'b', 'filled');

% Anomalie vere dei trasferimenti
% [~, ~, ~, ~, ~, theta_i_f] = car2par(ri_t_vect, vi_T_vect, mu); % Partenza
% [~, ~, ~, ~, ~, theta_f_t] = car2par(rf_t_vect, vf_A_vect, mu); % Arrivo

% Animazione della traiettoria
t = 0;

for k = 2:N

    [~, ~, ~, ~, ~, thetak_1] = car2par(r_1(:,k), v_1(:,k), mu);
    [~, ~, ~, ~, ~, thetak_2] = car2par(r_2(:,k), v_2(:,k), mu);
   
    % Aggiorna la traiettoria orbitale della Terra
    set(traj_1, 'XData', r_1(1,1:k), 'YData', r_1(2,1:k), 'ZData', r_1(3,1:k));
    % Aggiorna la posizione della Terra
    set(corpo_1, 'XData', r_1(1,k), 'YData', r_1(2,k), 'ZData', r_1(3,k)); 

    if((thetak_1 > thetai_1) && (thetaf_2 > thetak_2))
        t = t+1;
        % Aggiorna la traiettoria orbitale del satellite
        set(sat, 'XData', r_s(1,1:t), 'YData', r_s(2,1:t), 'ZData', r_s(3,1:t));
        % Aggiorno la posizione del satellite
        % Aggiorna le posizioni dei corpi
        set(satellite, 'XData', r_s(1,t), 'YData', r_s(2,t), 'ZData', r_s(3,t));
    end

    % Aggiorno la traiettoria e la posizione dell'asteroide
    set(traj_2, 'XData', r_2(1,1:k), 'YData', r_2(2,1:k), 'ZData', r_2(3,1:k));
    set(corpo_2, 'XData', r_2(1,k), 'YData', r_2(2,k), 'ZData', r_2(3,k));

    % Ruota lentamente la vista per una migliore prospettiva
    % camorbit(0.5,0);

    if(k == N)

        for j = N:(N+n-t)
            % Calcolo dell'anomalia vera successiva
            [~, theta2] = KepleroDiretto(a_1, e_1, theta_1(j), dt, mu);
            theta_1(j+1) = theta2;
            % Calcolo del vettore corrispondente
            [r_i, v_i] = par2car(a_1, e_1, i_1, OM_1, om_1, theta_1(j+1), mu);
            r_1(:, j+1) = r_i;
            v_1(:, j+1) = v_i;

            % Calcolo dell'anomalia vera successiva
            [~, theta2] = KepleroDiretto(a_2, e_2, theta_2(j), dt, mu);
            theta_2(j+1) = theta2;
            % Calcolo del vettore corrispondente
            [r_i, v_i] = par2car(a_2, e_2, i_2, OM_2, om_2, theta_2(j+1), mu);
            r_2(:, j+1) = r_i;
            v_2(:, j+1) = v_i;            
        end

        for y = N:(N+n-t)
            % Aggiorna la traiettoria orbitale del satellite
            set(sat, 'XData', r_s(1,1:y-(N-t)), 'YData', r_s(2,1:y-(N-t)), 'ZData', r_s(3,1:y-(N-t)));
            % Aggiorna le posizioni del satellite
            set(satellite, 'XData', r_s(1,y-(N-t)), 'YData', r_s(2,y-(N-t)), 'ZData', r_s(3,y-(N-t)));

            % Aggiorno la traiettoria e la posizione della Terra
            set(traj_1, 'XData', r_1(1,1:y), 'YData', r_1(2,1:y), 'ZData', r_1(3,1:y));
            set(corpo_1, 'XData', r_1(1,y), 'YData', r_1(2,y), 'ZData', r_1(3,y));

            % Aggiorno la traiettoria e la posizione dell'asteroide
            set(traj_2, 'XData', r_2(1,1:y), 'YData', r_2(2,1:y), 'ZData', r_2(3,1:y));
            set(corpo_2, 'XData', r_2(1,y), 'YData', r_2(2,y), 'ZData', r_2(3,y));

            pause(10e-8); % Tempo di pausa per il rendering fluido
        end

        % m = 0;
        % for j = (N+n-t):N_tot
        %     m = m+1;
        %     % Calcolo dell'anomalia vera successiva
        %     [~, theta2] = KepleroDiretto(a_1, e_1, theta_1(j), dt, mu);
        %     theta_1(j+1) = theta2;
        %     % Calcolo del vettore corrispondente
        %     [r_i, v_i] = par2car(a_1, e_1, i_1, OM_1, om_1, theta_1(j+1), mu);
        %     r_1(:, j+1) = r_i;
        %     v_1(:, j+1) = v_i;
        % 
        %     % Calcolo dell'anomalia vera successiva
        %     [~, theta2] = KepleroDiretto(a_2, e_2, theta_2(j), dt, mu);
        %     theta_2(j+1) = theta2;
        %     % Calcolo del vettore corrispondente
        %     [r_i, v_i] = par2car(a_2, e_2, i_2, OM_2, om_2, theta_2(j+1), mu);
        %     r_2(:, j+1) = r_i;
        %     v_2(:, j+1) = v_i;
        % 
        %     % Aggiorno la traiettoria e la posizione della Terra
        %     set(traj_1, 'XData', r_1(1,1:j), 'YData', r_1(2,1:j), 'ZData', r_1(3,1:j));
        %     set(corpo_1, 'XData', r_1(1,j), 'YData', r_1(2,j), 'ZData', r_1(3,j));
        % 
        %     % Aggiorno la traiettoria e la posizione dell'asteroide
        %     set(traj_2, 'XData', r_2(1,1:j), 'YData', r_2(2,1:j), 'ZData', r_2(3,1:j));
        %     set(corpo_2, 'XData', r_2(1,j), 'YData', r_2(2,j), 'ZData', r_2(3,j));
        % 
        %     pause(10e-7); % Tempo di pausa per il rendering fluido
        % end

    end
    
    pause(10e-8); % Tempo di pausa per il rendering fluido
end

fprintf('\n Intervalli temporali di discretizzazione totali : %d \n', N+n-t)
fprintf('\n Intervalli temporali di discretizzazione dalla partenza del satellite : %d \n', N_tot-t)
% fprintf('\n Intervalli temporali di discretizzazione del trasferimento : %d \n', N_tot-t-m)

% hold off;

end

