function plotorbit_3D_realSpeed(a, e, i, OM, om, th0, thf, dt, mu)

% Visualizzazione animata dell'orbita in 3D con una scia persistente con
% un'animazione dinamicamente coerente 
%
% ------------------------------------------------------------------
% Argomenti di Input
% a              Semi-asse maggiore                       [km]
% e              Eccentricità                             [-]
% i              Inclinazione                             [rad]
% OM             Ascensione retta del nodo ascendente     [rad]
% om             Anomalia del pericentro                  [rad]
% th0            Anomalia vera iniziale                   [rad]
% thf            Anomalia vera finale                     [rad]
% dt             Discretizzazione temporale               [s]
% mu             Costante gravitazionale                  [km^3/s^2]
%
% ------------------------------------------------------------------

% Tempo per arrivare all'anomalia iniziale dal perigeo
Ti = TOF(a, e, 0, th0, mu);
% Tempo orbitale dall'anomalia iniziale 
T = TOF(a, e, th0, thf, mu);
% Tempo complessivo dal perigeo
Tf = T + Ti;

% Vettore dei tempi
t_vect = Ti:dt:Tf;

N = length(t_vect);

% Inizializzazione del vettore delle anomalie vere
theta = zeros(1, N);
theta(1) = th0;

% Calcolo posizione ad istanti temporali regolari 
r = zeros(3, N);
[r_i, ~] = par2car(a, e, i, OM, om, th0, mu);
r(:, 1) = r_i;

for j = 1:N-1
    % Calcolo dell'anomalia vera successiva
    [~, theta2] = KepleroDiretto(a, e, theta(j), dt, mu);
    theta(j+1) = theta2;
    % Calcolo del vettore corrispondente
    [r_i, ~] = par2car(a, e, i, OM, om, theta(j+1), mu);
    r(:, j+1) = r_i;
end

% Adesso le posizioni delle posizioni sono equidistanziate temporalmente
% quindi non c'è l'esigenza di una pausa dinamica nell'animazione

% Creazione della figura
hold on
grid on
axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
title('3D Strategy Animation')
view(3)
rotate3d on

% Inizializzazione della traiettoria
traj = plot3(r(1,1), r(2,1), r(3,1), 'LineWidth', 1);

% Punto del satellite in movimento
satellite = scatter3(r(1,1), r(2,1), r(3,1), 25, 'r', 'filled');

% Animazione della traiettoria
for k = 2:N
    % Aggiorna la traiettoria
    set(traj, 'XData', r(1,1:k), 'YData', r(2,1:k), 'ZData', r(3,1:k));
    
    % Aggiorna la posizione del satellite
    set(satellite, 'XData', r(1,k), 'YData', r(2,k), 'ZData', r(3,k));
    
    % Ruota lentamente la vista per una migliore prospettiva
    % camorbit(0.5,0);
    
    pause(0.05); % Tempo di pausa per il rendering fluido
end

end


