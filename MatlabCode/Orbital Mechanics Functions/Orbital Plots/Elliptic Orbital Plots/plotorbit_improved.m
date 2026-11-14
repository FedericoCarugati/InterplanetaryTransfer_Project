function plotorbit_improved(a, e, i, OM, om, th0, thf, delta_th, mu)

% Visualizzazione animata dell'orbita in 3D con una scia persistente
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
% delta_th       Passo dell'anomalia vera                 [rad]
% mu             Costante gravitazionale                  [km^3/s^2]
%
% ------------------------------------------------------------------

% Controllo sugli input
if a <= 0 || e < 0 || e >= 1 || mu <= 0
    error('Parametri orbitali non validi.');
end

% Creazione del vettore delle anomalie vere
if thf < th0
    thf = thf + 2*pi;
end
theta = th0:delta_th:thf;
n = length(theta);

% Inizializzazione vettori per posizione
r = zeros(3, n);

% Calcolo delle posizioni orbitali
for j = 1:n
    r(:,j) = par2car(a, e, i, OM, om, theta(j), mu);
end

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
traj = plot3(r(1,1), r(2,1), r(3,1), 'LineWidth', 0.5);

% Punto del satellite in movimento
satellite = scatter3(r(1,1), r(2,1), r(3,1), 25, 'r', 'filled');

% Animazione della traiettoria
for k = 2:n
    % Aggiorna la traiettoria
    set(traj, 'XData', r(1,1:k), 'YData', r(2,1:k), 'ZData', r(3,1:k));
    
    % Aggiorna la posizione del satellite
    set(satellite, 'XData', r(1,k), 'YData', r(2,k), 'ZData', r(3,k));
    
    % Ruota lentamente la vista per una migliore prospettiva
    % camorbit(0.5,0);
    
    pause(0.002); % Tempo di pausa per il rendering fluido
end

% hold off;

end
