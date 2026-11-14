function delta_t = TOF(a, e, th1, th2, mu)

% Time of Flight
% 
% deltat = TOF(a, e, th1, th2, mu)
% 
% ------------------------------------------------------------------------
% Argomento di Input:
% a               [1x1]      semi-asse maggiore                [km]
% e               [1x1]      eccentricità                      [-]
% th1             [1x1]      anomalia vera iniziale            [rad]
% th2             [1x1]      anomalia vera finale              [rad]
% mu              [1x1]      costante gravitazional            [km^3/s^2]
% 
% ------------------------------------------------------------------------
% Argomento di Input:
% deltat          [1x1]      tempo di volo                     [s]
%
% ------------------------------------------------------------------------

% p = a*(1 - e^2);
% h = sqrt(p*mu);
% n = sqrt(h^3/mu^2);

n = sqrt(a^3/mu);
T = 2*sqrt(a^3/mu)*pi;

if(th1==th2)
    delta_t = 0;
else
    % Anomalie eccentriche
    E1 = 2*atan2(sqrt(1 - e) * sin(th1 / 2), sqrt(1 + e) * cos(th1 / 2));
    E2 = 2*atan2(sqrt(1 - e) * sin(th2 / 2), sqrt(1 + e) * cos(th2 / 2));
    
    deltat1 = n*(E1 - e*sin(E1));
    deltat2 = n*(E2 - e*sin(E2));
    
    % Tempo di volo
    if(th2 > th1)
        if(th2 > 2*pi)
            deltat2 = deltat2 + T;
        end
        if(th1 > 2*pi)
            deltat1 = deltat1 + T;
        end
        delta_t = deltat2 - deltat1;
    else
        if(th2 > 2*pi)
            deltat2 = deltat2 + T;
        end
        if(th1 > 2*pi)
            deltat1 = deltat1 + T;
        end
        delta_t = deltat2 - deltat1 + T;
    end
end

end

