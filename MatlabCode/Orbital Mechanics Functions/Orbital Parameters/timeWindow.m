function Deltat = timeWindow(a, e, theta1, Deltas_max, mu)

% Calcolo dell'intervallo temporale associato ad un intervallo spaziale nel
% moto orbitale di un corpo celeste
%
% -----------------------------------------------------------------------------
% Imput:
% a               [1x1]  semi asse maggiore                          [km]
% e               [1x1]  eccentricità                                [-]
% theta1          [1x1]  anomalia vera iniziale dell'intervallo      [rad]
% Deltas_max      [1x1]  Intervallo spaziale richiesto               [km]
% mu              [1x1]  parametro gravitazionale                    [km^3/s^2]
%
% ------------------------------------------------------------------------------
%Outputs:
% Deltat          [1x1]  Intervallo temporale associato              [s]
%
% ------------------------------------------------------------------------------

% La funzione timeWindow serve a calcolare l'intervallo temporale
% necessario per avere un massimo intervallo spaziale associato,
% considerando il moto orbitale di un corpo celeste

% Nello specifico questa funzione è stata adattata in modo da identificare
% il corretto intervallo temporale da sfruttare per il calcolo del
% trasferimento interplanetario ottimale dello Scenario 02

dt = linspace(0.1, 5, 1e2);
n = length(dt);
Deltas_vect = zeros(n,1);
theta2_vect = zeros(n,1);

% Si parametrizza l'ellisse
% x = a*cos(theta)
% y = b*sin(theta)
% Dunque l'arco sarà l'integrale rispetto a theta di sqrt(x^2 + y^2)

% Calcolo del semiasse minore dell'orbita dell'asteroide
b = a * sqrt(1 - e^2);

max = 0;

for(j = 1:n)
    [~, theta2] = KepleroDiretto(a, e, theta1, dt(j), mu);
    theta2_vect(j) = theta2;

    % Conversione da anomalia vera a anomalia eccentrica (E)
    E1 = 2 * atan2(sqrt(1 - e) * sin(theta1 / 2), sqrt(1 + e) * cos(theta1 / 2));
    E2 = 2 * atan2(sqrt(1 - e) * sin(theta2 / 2), sqrt(1 + e) * cos(theta2 / 2));
    
    % Correzione per dominio dell'angolo
    if (E1 < 0)
        E1 = E1 + 2*pi;
    end
    if (E2 < 0)
        E2 = E2 + 2*pi;
    end
    
    % Funzione integranda
    integranda = @(E) sqrt(a^2*sin(E).^2 + b^2*cos(E).^2);
    
    % Calcolo dell'integrale numerico tra E1 ed E2
    Deltas_vect(j) = abs(integral(integranda, E1, E2));  % Lunghezza in km

    if(Deltas_vect(j) > Deltas_max && max == 0)
        max = j;
    end
end

Deltat = dt(max-1);

end

