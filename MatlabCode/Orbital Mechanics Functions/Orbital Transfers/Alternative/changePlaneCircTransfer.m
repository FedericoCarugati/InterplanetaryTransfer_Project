function [deltaV_tot, deltaV_circ_1, deltaV_trasf, deltaV_circ_2, deltaV_ecc, deltaV_plane, ...
    deltaT_circ_1, deltaT_circ_2, deltaT_t, deltaT, a_t, e_t, om_t, theta_circ_1_ell, ...
    theta_circ_2, theta_circ_2_ell] = changePlaneCircTransfer(a_i, e_i, om_i, a_f, e_f, ...
    om_f, mu, OM_i, i_i, OM_f, i_f, theta_t)

% Trasferimento orbitale di cambio di piano con circolarizzazione
% (Approccio "semi-biellittico")
% 
% [deltaV_tot, deltaV_circ_1, deltaV_trasf, deltaV_circ_2, deltaV_ecc, ...
%     deltaT_circ_1, deltaT_circ_2, deltaT_t, deltaT, a_t, e_t, om_t, theta_circ_1_ell, ...
%     theta_circ_2, theta_circ_2_ell] = changePlaneCircTransfer(a_i, e_i, om_i, a_f, e_f, ...
%     om_f, mu, OM_i, i_i, OM_f, i_f, theta_t)
% 
% ----------------------------------------------------------------------------
% Argomenti di input:
% a_i               [1x1]     semi-asse maggiore iniziale               [km]
% e_i               [1x1]     eccentricità iniziale                     [-]
% om_i              [1x1]     anomalia di pericentro orbita iniziale    [rad]
% a_f               [1x1]     semi-asse maggiore finale                 [km]
% e_f               [1x1]     eccentricità finale                       [-]
% om_f              [1x1]     anomalia di pericentro orbita finale      [rad]
% mu                [1x1]     costante gravitazionale                   [km^3/s^2]
% OM_i              [1x1]     RAAN iniziale                             [rad]
% i_i               [1x1]     inclinazione iniziale dell'orbita         [rad]
% OM_f              [1x1]     RAAN finale                               [rad]
% i_f               [1x1]     inclinazione finale dell'orbita           [rad]
% theta_t           [1x1]     anomalia vera prima manovra biellittica   [rad]
% 
% ----------------------------------------------------------------------------
% Argomenti di output:
% DeltaV_tot        [1x1]     costo complessivo della manovra           [km/s]
% DeltaV_circ_1     [1x1]     costo della circolarizzazione 1           [km/s]
% deltaV_trasf      [1x1]     costo del trasferimento ellittico         [km/s]
% DeltaV_circ_2     [1x1]     costo della circolarizzazione 2           [km/s]
% Deltat_ecc        [1x1]     costo dell'ellitticizzazione              [km/s]
% DeltaV_plane      [1x1]     impulso del cambio di piano               [km/s]
% Deltat_circ_1     [1x1]     tempo sull'orbita circolare 1             [s]
% Deltat_circ_2     [1x1]     tempo sull'orbita circolare 2             [s]
% DeltaT_t          [1x1]     tempo sull'ellisse di trasferimento       [s]
% DeltaT            [1x1]     tempo totale del trasferimento            [s]
% a_t               [1x1]     semi-asse orbita di trasferimento         [km]
% e_t               [1x1]     eccentricità orbita di trasferimento      [-]
% om_t              [1x1]     anomalia del pericentro di trasferimento  [rad]
% theta_circ_1_ell  [1x1]     anomalia vera orbita circolare 1 nel trasf[rad]
% theta_circ_2      [1x1]     anomalia vera orbita circolare 2 finale   [rad]
% theta_circ_2_ell  [1x1]     anomalia vera orbita circolare 2 nel trasf[rad]
%
% ----------------------------------------------------------------------------

% Questa manovra si basa su una iniziale circolarizzazione dell'orbita
% iniziale all'apogeo per poi iniziare la manovra di cambio di piano nel punto
% di intersezione tra i due piani orbitali
% Nel punto di cambio di piano si opera anche una manovra in cui si rende
% l'orbita ellittica con parametri orbitali stabiliti dalle dimensioni delle
% orbite ellittiche iniziale e finale
% L'ultimo step della manovra è quello di rircolarizzare l'orbita con un
% raggio dato dalla distanza all'apogeo dell'orbita finale per arrivare
% all'apogeo di questa e inserirsi in quella traiettoria

% Parametri orbitali dell'orbita inziale
% rp_i = a_i*(1 - e_i); % Pericentro dell'orbita iniziale
ra_i = a_i*(1 + e_i); % Apocentro dell'orbita iniziale

% Parametri orbitali dell' orbita finale
% rp_f = a_f*(1 - e_f); % Pericentro dell'orbitale finale
ra_f = a_f*(1 + e_f); % Apocentro dell'orbitale finale

if(ra_i > ra_f)
    ra_t = ra_i; % apogeo orbita di trasferimento
    rp_t = ra_f; % perigeo orbita di trasferimento
else
    ra_t = ra_f; % apogeo orbita di trasferimento
    rp_t = ra_i; % perigeo orbita di trasferimento
end

% Parametri dell'orbita ellittica di trasferimento
e_t = (ra_t - rp_t)/(ra_t + rp_t); % Eccentricità dell'orbita di trasferimento
a_t = (ra_t + rp_t)/2; % Semiasse maggiore dell'orbita di trasferimento
p_t = a_t*(1 - e_t^2); % Semilato retto dell'orbita iniziale
va_t = sqrt(mu/p_t)*(1 - e_t); % Velocità all'apogeo dell'orbita di trasferimento
% vp_t = sqrt(mu/p_t)*(1 + e_t); % Velocità al perigeo dell'orbita di trasferimento
vp_t = sqrt(mu*(2/rp_t - 1/a_t));

p_i = a_i*(1 - e_i^2); % Semilato retto dell'orbita iniziale
p_f = a_f*(1 - e_f^2); % Semilato retto dell'orbita iniziale

% Parametri delle orbite circolari
r_circ_1 = ra_i;
r_circ_2 = ra_f;

% Calcolo dell'impulso per circolarizzare l'orbita iniziale
va_i = sqrt(mu/p_i)*(1 - e_i); % Velocità dell'orbita iniziale all'apogeo
v_circ_1 = sqrt(mu/r_circ_1); % Velocità prima orbita circolare

deltaV_circ_1 = abs(va_i - v_circ_1); % Costo manovra di circolarizzazione

% Questa manovra si opera all'anomalia vera theta_i = pi rad dell'orbita
% iniziale e corrisponde all'anomalia vera theta_circ_1 = (om_i + pi) rad
% dell'orbita circolare

% Quindi dall'apogeo dell'orbita ellittica iniziale si passa ad un'orbita
% circolare che prosegue il proprio cammino fino al punto in cui si passa
% all'orbita di trasferimento ellittica

% A questo punto si distinguono due casi a seconda delle dimensioni
% iniziale e finale delle orbite

% Infatti il trasferimento di cambio di piano è ottimale quando la velocità
% orbitale è minore:

% r_circ_1 = ra_i e r_circ_2 = ra_f, dunque la manovra sarà svolta nel
% punto più lontano dall'attrattore e, inoltre, è da notare che:
% Se r_circ = r_a => v_circ > v_a
% Se r_circ = r_p => v_p > v_circ

if(ra_f > ra_i)
    % La velocità minore si avrà nel punto di apogeo dell'orbita di
    % trasferimento (nel punto della seconda circolarizzazione) e la
    % velocità sull'ellisse in questo punto è minore della corrispondente
    % velocità sull'orbita circolare alla stessa distanza dall'attratore

    % La manovrà si svolge sull'orbita ellittica di trasferimento quindi si
    % calcolano i costi di trasferimento da operare per raggiungere questa
    % orbita a partire dalla prima orbita circolare

    % Ellitticizzazione dell'orbita
    deltaV_trasf = abs(v_circ_1 - vp_t);

    % ----------------------------------------------------------------------------
    % Anomalia vera della manovra sull'orbita circolare:

    % Se theta_t + om_i > 2*pi:
    % theta_circ_1_ell = pi + (om_i - (2*pi - theta_t)) 

    % Se 2*pi > theta_t + om_i && theta_t > pi:
    % theta_circ_1_ell = pi + (2*pi - (om_i + theta_t))

    % Se 2*pi > theta_t + om_i && pi > theta_t:
    % theta_circ_1_ell = pi + (om_i + theta_t) 
    % ----------------------------------------------------------------------------

    if(theta_t + om_i >= 2*pi)
        theta_circ_1_ell = pi + (om_i - (2*pi - theta_t));
        fprintf('theta_t + om_i >= 2*pi \n')
    else
        if(theta_t >= pi)
            theta_circ_1_ell = pi + (2*pi - (om_i + theta_t));
            fprintf('theta_t >= pi \n')
        else
            theta_circ_1_ell = (om_i + theta_t) - pi;
            fprintf('pi > theta_t \n')
        end
    end

    % Anomalia vera della manovra sull'orbita ellittica di trasferimento:
    % theta_t_i = 0 (perigeo dell'orbita di trasferimento)

    % A questo punto si procede sulla nuova orbita ellittica fino
    % all'apogeo
    % Dunque si svolgerà la manovra di cambio di piano seguita dalla
    % manovra di circolarizzazione dell'orbita

    % Manovra di cambio piano
    [deltaV_plane, om_t, theta_plane, ~] = changeOrbitalPlane(a_t, e_t, i_i, OM_i, theta_circ_1_ell, i_f, OM_f, mu);
    % Verifica del cambio di piano
    if(theta_plane ~= pi)
        fprintf('Problema con il punto di cambio piano (anomalia) : %d \n', theta_plane)
    end
    % Se theta_plane == pi allora la manovra di cambio piano è stata
    % correttamente svolta al perigeo delle due orbite ellittiche (la
    % seconda viene subito trasformata in un'orbita circolare) 

    v_circ_2 = sqrt(mu/r_circ_2); % Velocità seconda orbita circolare

    % Manovra di circolarizzazione all'apogeo
    deltaV_circ_2 = abs(va_t - v_circ_2);

    % Anomalia vera nella seconda circolarizzazione rispetto all'orbita circolare
    % theta_circ_2 = om_t - pi;
    theta_circ_2 = om_t + pi;

else
    % Questa parte dell'if è completamente da rivedere !!!!

    % In questo caso la velocità minore si avrà nel punto di apogeo dell'orbita di
    % trasferimento (nel punto della prima circolarizzazione) e la
    % velocità sull'ellisse in questo punto è minore della corrispondente
    % velocità sull'orbita circolare alla stessa distanza dall'attratore

    % Dunque si svolgerà la manovra di ellitticizzazione dell'orbita
    % seguita dalla manovra di cambio di piano

    % Ellitticizzazione dell'orbita
    deltaV_trasf = abs(v_circ_1 - vp_t);

    % Anomalia vera ellitticizzazione per l'orbita circolare
    if(theta_t + om_i >= 2*pi)
        theta_circ_1_ell = (om_i - (2*pi - theta_t));
    else
        if(theta_t >= pi)
            theta_circ_1_ell = (2*pi - (om_i + theta_t));
        else
            theta_circ_1_ell = (om_i + theta_t);
        end
    end

    % Manovra di cambio piano
    [deltaV_plane, om_t, theta_plane, ~] = changeOrbitalPlane(a_t, e_t, i_i, OM_i, om_i, i_f, OM_f, mu);
    if(om_t ~= theta_t || theta_plane ~= pi)
        fprintf('Problema con il cambio piano')
        % fprintf('%d; %d;' & om_t, theta_plane)
    end

    % Infine si passa alla manovra di circolarizzazione dell'orbita
    % ellittica

    v_circ_2 = sqrt(mu/r_circ_2); % Velocità seconda orbita circolare

    % Manovra di circolarizzazione all'apogeo
    deltaV_circ_2 = abs(va_t - v_circ_2);

    % Anomalia vera nella seconda circolarizzazione rispetto all'orbita circolare
    % theta_circ_2 = theta_circ_1_ell - pi;
    theta_circ_2 = om_t - pi;
end

% Il satellite viaggerà sull'orbita circolare fino ad arrivare al punto di
% apogeo dell'orbita finale in cui opererà l'eccentricizzazione

% Il punto per operare questa manovra viene preso in modo tale da ottenere
% un'orbita ellittica finale che abbia già l'anomalia di pericentro finale
% desiderata

% Il satellite percorre l'orbita circolare per un angolo calcolato rispetto
% al piano equatoriale di (pi + om_f) gradi e arriva all'anomalia vera pi
% dell'orbita finale (rispetto all'anomalia di pericentro dell'orbita
% finale)

% Anomalia vera percorso dal satellite sulla seconda orbita circolare 
% rispetto al piano equatoriale terrestre
theta_circ_2_ell = pi + om_f; 

va_f = sqrt(mu/p_f)*(1 - e_f); % Velocità dell'orbita finale all'apogeo
deltaV_ecc = abs(va_f - v_circ_2); % Costo manovra di eccentricizzazione

% Calcolo dei costi di trasferimento
deltaV_tot = deltaV_circ_1 + deltaV_circ_2 + deltaV_plane + ...
    deltaV_trasf + deltaV_ecc;

% Calcolo dei costi della manovra
deltaT_circ_1 = (theta_circ_1_ell - pi)/v_circ_1;
deltaT_circ_2 = (om_f + pi - theta_circ_2)/v_circ_2;
deltaT_t = pi*sqrt(a_t^3/mu);
deltaT = deltaT_circ_1 + deltaT_circ_2 + deltaT_t;

end

