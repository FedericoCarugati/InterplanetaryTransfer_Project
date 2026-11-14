function [delta_theta, theta2] = KepleroDiretto(a, e, theta1, deltaT, mu)

% Questa funzione serve per il problema diretto della legge oraria per
% orbita ellittica e offre la differenza tra le anomalie vere data
% l'anomalia iniziale e il tempo del moto

% KepleroInvero
% 
% [delta_theta, theta2] = KepleroInverso(a, e, theta1, deltaT, mu)
% 
% ------------------------------------------------------------------------
% Argomento di Input:
% a               [1x1]      semi-asse maggiore                [km]
% e               [1x1]      eccentricità                      [-]
% theta1          [1x1]      anomalia vera iniziale            [rad]
% deltaT          [1x1]      tempo tra le anomalie             [s]
% mu              [1x1]      costante gravitazional            [km^3/s^2]
% 
% ------------------------------------------------------------------------
% Argomento di Input:
% delta_theta     [1x1]      Differenza tra le anomalie vere   [rad]
% theta2          [1x1]      Anomalia vera finale              [rad]
%
% ------------------------------------------------------------------------

% Per il calcolo corretto della seconda anomalia serve ricavare il tempo di
% volo per raggiungere la prima anomalia dal punto di perigeo dell'robita
% Ottenuto questo valore, sarà possibile calcolare il tempo complessivo per
% giungere alla seconda anomalia vera dal punto di pericentro e, dunque,
% calcolare il suo valore in radianti

t1 = TOF(a, e, 0, theta1, mu);

% Il problema diretto richiede la soluzione di un'equazione implicita dunque
% si implementerà anche la funzione di newton per la risoluzione di sistemi 
% non lineari

% Costante associata al problema 
n = sqrt(mu/a^3);

% Anomalia media 
M = n*(deltaT + t1);

% Il calcolo dell'anomalia eccentrica richiede la risoluzione di
% un'equazione implicita (risoluzione tramite Metodo di Newton)

% Definizione delle opzioni per fsolve
options = optimoptions('fsolve', 'Display', 'off');

% Funzione di Newton
fun = @(x) x - e*sin(x) - M;
x0 = 1;
E = fsolve(fun, x0, options);

% Anomalia eccentrica
if(0 >= E)
    E = 2*pi + E;
end

% Anomalia vera
theta_next = 2*atan2(sqrt(1 + e)*sin(E/2), sqrt(1 - e)*cos(E/2));
if(theta_next >= 0)
    theta2 = theta_next;
else
    theta2 = 2*pi + theta_next;
end

% Differenza tra le anomalie
delta_theta = theta2 - theta1;

end

