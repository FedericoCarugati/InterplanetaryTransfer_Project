function [DeltaV_tot, DeltaV1, DeltaV2, DeltaV_circ_1, DeltaV_circ_2, ...
    DeltaV_plane, DeltaV_ecc, Deltat1, Deltat2, DeltaT, a_1, e_1, om_1, ...
    om_2, a_2, e_2, theta_plane1, gamma, delta, v] = biellipticTransferChangePlane(...
    a_i, e_i, om_i, a_f, e_f, om_f, ra_t, mu, OM_i, i_i, OM_f, i_f, theta_t)

% Trasferimento orbitale bi-ellittico con cambio di piano
% 
% [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2, DeltaT, theta,...
%     a_1, e_1, om_1, theta_plane1, theta_plane2, om_f] = ...
%     biellipticTransferChangePlane(r_vect_i, v_vect_i, a_i, e_i, a_f,...
%     e_f, ra_t, mu, OM_f, i_f)
% 
% ----------------------------------------------------------------------------
% Argomenti di input:
% a_i               [1x1]     semi-asse maggiore iniziale               [km]
% e_i               [1x1]     eccentricità iniziale                     [-]
% om_i              [1x1]     anomalia di pericentro orbita iniziale    [rad]
% a_f               [1x1]     semi-asse maggiore finale                 [km]
% e_f               [1x1]     eccentricità finale                       [-]
% om_f              [1x1]     anomalia di pericentro orbita finale      [rad]
% ra_t              [1x1]     apocentro dell'orbita di trasferimento    [km]
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
% DeltaV1           [1x1]     impulso della prima manovra               [km/s]
% DeltaV2           [1x1]     impulso della seconda manovra             [km/s]
% DeltaV_circ_1     [1x1]     costo della circolarizzazione             [km/s]
% DeltaV_circ_2     [1x1]     costo della circolarizzazione             [km/s]
% DeltaV_plane      [1x1]     costo del cambio di piano orbitale        [km/s]
% DeltaV_ecc        [1x1]     costo dell'eccentricizzazione finale      [km/s]
% Deltat1           [1x1]     tempo della prima manovra                 [s]
% Deltat2           [1x1]     tempo della seconda manovra               [s]
% DeltatT           [1x1]     tempo totale del trasferimento            [s]
% om_f              [1x1]     anomalia del pericentro finale            [km]
% a_1               [1x1]     semi-asse orbita di trasferimento         [km]
% e_1               [1x1]     eccentricità orbita di trasferimento      [-]
% om_1              [1x1]     anomalia del pericentro orbita di trasf   [rad]
% om_2              [1x1]     anomalia del pericentro orbita di trasf   [rad]
% a_2               [1x1]     semi-asse orbita di trasferimento         [km]
% e_2               [1x1]     eccentricità orbita di trasferimento      [-]
% theta_plane1      [1x1]     anomalia vera sul primo trasferimento     [rad]
% gamma             [1x1]     anomalia di pericentro orbite di trasf    [rad]
% delta             [1x1]     anomalia vera orbita circolare finale     [rad]
% v                 [1x1]     scelta strategica ottimale                [-]
%
% ----------------------------------------------------------------------------

% Questa manovra si basa su una iniziale circolarizzazione dell'orbita
% iniziale all'apogeo per poi iniziare la manovra biellittica nel punto di
% intersezione con l'orbita finale. 

% Parametri orbitali dell'orbita inziale
% rp_i = a_i*(1 - e_i); % Pericentro dell'orbita iniziale
ra_i = a_i*(1 + e_i); % Apocentro dell'orbita iniziale

% Parametri orbitali dell' orbita finale
% rp_f = a_f*(1 - e_f); % Pericentro dell'orbitale finale
ra_f = a_f*(1 + e_f); % Apocentro dell'orbitale finale

% Prametri orbitali delle orbite di trasferimento

if(ra_t > ra_f && ra_f > ra_i || ra_t > ra_i && ra_i > ra_f)
    rp_1 = ra_i;  % Pericentro della prima orbita di trasferimento
    rp_2 = ra_f; % Pericentro della seconda orbita di trasferimento

    ra_1 = ra_t; % Apocentro della prima orbita di trasferimento
    ra_2 = ra_t; % Apocentro della seconda orbita di trasferimento
end

if(ra_f > ra_t && ra_t > ra_i)
    rp_1 = ra_i;  % Pericentro della prima orbita di trasferimento
    rp_2 = ra_t; % Pericentro della seconda orbita di trasferimento

    ra_1 = ra_t; % Apocentro della prima orbita di trasferimento
    ra_2 = ra_f; % Apocentro della seconda orbita di trasferimento
end

if(ra_f > ra_t && ra_i > ra_t || ra_i > ra_f && ra_f > ra_t)
    rp_1 = ra_t;  % Pericentro della prima orbita di trasferimento
    rp_2 = ra_t; % Pericentro della seconda orbita di trasferimento

    ra_1 = ra_i; % Apocentro della prima orbita di trasferimento
    ra_2 = ra_f; % Apocentro della seconda orbita di trasferimento
end

if(ra_i > ra_t && ra_t > ra_f)
    rp_1 = ra_t;  % Pericentro della prima orbita di trasferimento
    rp_2 = ra_f; % Pericentro della seconda orbita di trasferimento

    ra_1 = ra_i; % Apocentro della prima orbita di trasferimento
    ra_2 = ra_t; % Apocentro della seconda orbita di trasferimento
end

a_1 = (ra_1 + rp_1)/2; % semi-asse maggiore della prima orbita di trasferimento
a_2 = (ra_2 + rp_2)/2; % semi-asse maggiore della seconda orbita di trasferimento

p_i = a_i*(1 - e_i^2); % Semilato retto dell'orbita iniziale
p_f = a_f*(1 - e_f^2); % Semilato retto dell'orbita iniziale

e_1 = (ra_1 - rp_1)/(ra_1 + rp_1); % eccentricità della prima orbita di trasferimento
e_2 = (ra_2 - rp_2)/(ra_2 + rp_2); % eccentricità della seconda orbita di trasferimento

p_1 = a_1*(1 - e_1^2); % Semilato retto della seconda orbita di trasferimento
p_2 = a_2*(1 - e_2^2); % Semilato retto della seconda orbita di trasferimento

% Parametri delle orbite circolari
r_circ_1 = ra_i;
r_circ_2 = ra_f;

% Calcolo dell'impulso per circolarizzare l'orbita iniziale
va_i = sqrt(mu/p_i)*(1 - e_i); % Velocità dell'orbita iniziale all'apogeo
v_circ_1 = sqrt(mu/r_circ_1); % Velocità prima orbita circolare

DeltaV_circ_1 = abs(va_i - v_circ_1); % Costo manovra di circolarizzazione

% Dall'apogeo dell'orbita ellittica iniziale si passa ad un'orbita
% circolare che prosegue il proprio cammino fino al punto di inizio di
% manovra biellittica

% A questo punto bisogna calcolare il punto in cui operare la manovra nel
% punto di intersezioni dei piani 

% Anomalia di inizio della manovra biellittica rispetto al nodo ascendente
if(ra_t > ra_f && ra_t > ra_i)
    gamma = om_i - (pi - theta_t);
    if(0 > gamma)
        gamma = 2*pi + gamma;
    end
    if(gamma > 2*pi)
        gamma = gamma - 2*pi;
    end
else
    gamma = om_i - (2*pi - theta_t);
    if(0 > gamma)
        gamma = 2*pi + gamma;
    end
end
% Angolo a cui avviene la prima manovra biellittica sull'orbita circolare
% rispetto al piano equatoriale terrestre (anomalia di pericentro per le
% orbite ellittiche del trasferimento biellittico)

if(ra_t > ra_f && ra_f > ra_i || ra_t > ra_i && ra_i > ra_f)
    v = 1;
    % Calcolo della velocità desiderata nel primo punto della manovra
    % biellittica per spostarsi sulla prima orbita ellittica
    vp_1 = sqrt(mu*(2/rp_1 - 1/a_1));
    
    % Costo dell'inizio di manovra biellittica
    DeltaV1 = abs(vp_1 - v_circ_1);

    % A questo punto si giunge sulla prima orbita di trasferimento su cui si
    % andrà ad operare un cambio di piano all'apogeo (o perigeo), dunque a pi gradi
    % rispetto alla prima manovra del trasferimento biellittico
    
    va_1 = sqrt(mu*(2/ra_1 - 1/a_1)); % Velocità all'apogeo della prima orbita di trasferimento
    va_2 = sqrt(mu*(2/ra_2 - 1/a_2)); % Velocità all'apogeo della seconda orbita di trasferimento

    % Il costo della manovra di cambio piano è influenzata dal modulo della
    % velocità dunque si sceglie di operare il cambio di piano nel punto
    % con la velocità minore (apogeo dell'orbita con maggiore eccentricità, 
    % ovvero la seconda)

    if(va_1 > va_2)
        % Costo della seconda manovra del trasferimento biellittico
        DeltaV2 = abs((va_1 - va_2));
        
        % Manovra di cambio di piano nel punto di apogeo delle orbite di 
        % trasferimento della manovra biellittica
        om_1 = gamma;
        [DeltaV_plane, om_2, theta_plane1, ~] = ...
            changeOrbitalPlane(a_2, e_2, i_i, OM_i, om_1, i_f, OM_f, mu);
    else     
        % Manovra di cambio di piano nel punto di apogeo delle orbite di 
        % trasferimento della manovra biellittica
        om_1 = gamma;
        [DeltaV_plane, om_2, theta_plane1, ~] = ...
            changeOrbitalPlane(a_1, e_1, i_i, OM_i, om_1, i_f, OM_f, mu);

        % Costo della seconda manovra del trasferimento biellittico
        DeltaV2 = abs((va_1 - va_2));
    end

    % A questo punto si arriva al perigeo (o apogeo) della seconda orbita di trasferimento 
    % biellittico (che sarà uguale all'apogeo dell'orbita finale) e in questo
    % punto si vuole operare una nuova circolarizzazione dell'orbita per
    % garantire il cambio finale nel punto di apogeo dell'orbita finale
    
    vp_2 = sqrt(mu/p_2)*(1 + e_2); % Velocità della seconda orbita di trasferimento al perigeo
    v_circ_2 = sqrt(mu/r_circ_2); % Velocità seconda orbita circolare 
    
    DeltaV_circ_2 = abs(vp_2 - v_circ_2); % Costo manovra di circolarizzazione
    
    % Il satellite viaggerà sull'orbita circolare fino ad arrivare al punto di
    % apogeo dell'orbita finale in cui opererà l'eccentricizzazione
    
    % Il punto per operare questa manovra viene preso in modo tale da ottenere
    % un'orbita ellittica finale che abbia già l'anomalia di pericentro finale
    % desiderata
    
    % Il satellite percorre l'orbita circolare per un angolo calcolato rispetto
    % al piano equatoriale di (pi + om_f) gradi e arriva all'anomalia vera pi
    % dell'orbita finale (rispetto all'anomalia di pericentro dell'orbita
    % finale)
end

if(ra_f > ra_t && ra_t > ra_i || ra_f > ra_t && ra_i > ra_t)
    v = 2;
    % Calcolo della velocità desiderata nel primo punto della manovra
    % biellittica per spostarsi sulla prima orbita ellittica
    vp_1 = sqrt(mu*(2/rp_1 - 1/a_1));
    
    % Costo dell'inizio di manovra biellittica
    DeltaV1 = abs(vp_1 - v_circ_1);

    % A questo punto si giunge sulla prima orbita di trasferimento su cui si
    % andrà ad operare un cambio di piano all'apogeo (o perigeo), dunque a pi gradi
    % rispetto alla prima manovra del trasferimento biellittico
    
    va_1 = sqrt(mu*(2/ra_1 - 1/a_1)); % Velocità all'apogeo della prima orbita di trasferimento
    vp_2 = sqrt(mu*(2/rp_2 - 1/a_2)); % Velocità al perigeo della seconda orbita di trasferimento
    
    % Costo della seconda manovra del trasferimento biellittico
    DeltaV2 = abs((va_1 - vp_2));

    % La manovra di cambio piano viene svolta con velocità minima quindi si 
    % svolge prima della circolarizzazione perché la velocità dell'orbita 
    % ellittica all'apogeo sarà minore di quella circolare alla stessa 
    % distanza dall'attratore
    
    % Manovra di cambio di piano nel punto di apogeo delle orbite di 
    % trasferimento della manovra biellittica
    om_1 = gamma;
    [DeltaV_plane, om_2, theta_plane1, ~] = ...
        changeOrbitalPlane(a_2, e_2, i_i, OM_i, om_1, i_f, OM_f, mu);
    
    va_2 = sqrt(mu/p_2)*(1 - e_2); % Velocità della seconda orbita di trasferimento all'apogeo
    v_circ_2 = sqrt(mu/r_circ_2); % Velocità seconda orbita circolare 
    
    DeltaV_circ_2 = abs(va_2 - v_circ_2); % Costo manovra di circolarizzazione
    
    % Il satellite viaggerà sull'orbita circolare fino ad arrivare al punto di
    % apogeo dell'orbita finale in cui opererà l'eccentricizzazione
    
    % Il punto per operare questa manovra viene preso in modo tale da ottenere
    % un'orbita ellittica finale che abbia già l'anomalia di pericentro finale
    % desiderata
    
    % Il satellite percorre l'orbita circolare per un angolo calcolato rispetto
    % al piano equatoriale di (pi + om_f) radianti e arriva all'anomalia vera pi
    % dell'orbita finale (rispetto all'anomalia di pericentro dell'orbita
    % finale)
end

if(ra_i > ra_f && ra_f > ra_t || ra_i > ra_t && ra_t > ra_f)
    v = 3;
    % In questo caso il cambio piano è la prima manovra ad essere operata e
    % viene svolta appena dopo l'ellitticizzazione

    % Calcolo della velocità desiderata nel primo punto della manovra
    % biellittica per spostarsi sulla prima orbita ellittica
    va_1 = sqrt(mu*(2/ra_1 - 1/a_1));
    
    % Costo dell'inizio di manovra biellittica
    DeltaV1 = abs(va_1 - v_circ_1);

    % Manovra di cambio di piano nel punto di apogeo delle orbite di 
    % trasferimento della manovra biellittica
    om_1 = gamma;
    [DeltaV_plane, om_2, theta_plane1, ~] = ...
        changeOrbitalPlane(a_1, e_1, i_i, OM_i, om_1, i_f, OM_f, mu);
    
    vp_1 = sqrt(mu*(2/rp_1 - 1/a_1)); % Velocità al perigeo della prima orbita di trasferimento
    vp_2 = sqrt(mu*(2/rp_2 - 1/a_2)); % Velocità al perigeo della seconda orbita di trasferimento
    
    % Costo della seconda manovra del trasferimento biellittico
    DeltaV2 = abs((vp_1 - vp_2));

    % Dunque si svolge l'ultima manovra impulsiva per circolarizzare
    % l'orbita
    
    va_2 = sqrt(mu/p_2)*(1 - e_2); % Velocità della seconda orbita di trasferimento all'apogeo
    v_circ_2 = sqrt(mu/r_circ_2); % Velocità seconda orbita circolare 
    
    DeltaV_circ_2 = abs(va_2 - v_circ_2); % Costo manovra di circolarizzazione
    
    % Il satellite viaggerà sull'orbita circolare fino ad arrivare al punto di
    % apogeo dell'orbita finale in cui opererà l'eccentricizzazione
    
    % Il punto per operare questa manovra viene preso in modo tale da ottenere
    % un'orbita ellittica finale che abbia già l'anomalia di pericentro finale
    % desiderata
    
    % Il satellite percorre l'orbita circolare per un angolo calcolato rispetto
    % al piano equatoriale di (pi + om_f) gradi e arriva all'anomalia vera pi
    % dell'orbita finale (rispetto all'anomalia di pericentro dell'orbita
    % finale)
end

% Angolo percorso dal satellite sulla seconda orbita circolare rispetto al
% piano equatoriale terrestre
delta = pi + om_f; 

% Calcolo dell'impulso per eccentricizzare la seconda orbita circolare
va_f = sqrt(mu/p_f)*(1 - e_f); % Velocità dell'orbita finale all'apogeo
v_circ_2 = sqrt(mu/r_circ_2); % Velocità della seconda orbita circolare

DeltaV_ecc = abs(va_f - v_circ_2); % Costo manovra di eccentricizzazione

% Calcolo dei costi di trasferimento
DeltaV_tot = DeltaV1 + DeltaV2 + DeltaV_circ_1 + DeltaV_circ_2 + ...
    DeltaV_plane + DeltaV_ecc;

% Calcolo dei tempi di trasferimento (delta_T)
Deltat_circ_1 = (abs(pi - theta_t)*ra_i)/(v_circ_1);
Deltat_circ_2 = (abs(delta)*ra_f)/(v_circ_2);
Deltat1 = pi*sqrt(a_1^3/mu);
Deltat2 = pi*sqrt(a_2^3/mu);

DeltaT = Deltat1 + Deltat2 + Deltat_circ_1 + Deltat_circ_2;

end

