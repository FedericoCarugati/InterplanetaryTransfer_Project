function [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2, DeltaT] = ...
    biellipticTransfer(a_i, e_i, a_f, e_f, ra_t, mu)

% Trasferimento orbitale biellittico
% 
% [DeltaV1, DeltaV2, DeltaV3, Deltat1, Deltat2, DeltaT] = ...
% biellipticTransfer(a_i, e_i, a_f, e_f, ra_f, mu)
% 
% ----------------------------------------------------------------------------
% Argomenti di input:
% a_i               [1x1]     semi-asse maggiore iniziale               [km]
% e_i               [1x1]     eccentricità iniziale                     [-]
% a_f               [1x1]     semi-asse maggiore finale                 [km]
% e_f               [1x1]     eccentricità finale                       [-]
% ra_t              [1x1]     apocentro dell'orbita di trasferimento    [km]
% mu                [1x1]     costante gravitazionale                   [km^3/s^2]
% 
% ----------------------------------------------------------------------------
% Argomenti di output:
% DeltaV1           [1x1]     impulso della prima manovra               [km/s]
% DeltaV2           [1x1]     impulso della seconda manovra             [km/s]
% DeltaV3           [1x1]     impulso della terza manovra               [km/s]
% Deltat1           [1x1]     tempo della prima manovra                 [s]
% Deltat2           [1x1]     tempo della seconda manovra               [s]
% DeltatT           [1x1]     tempo totale del trasferimento            [s]
%

% Parametri orbitali dell'orbita inziale
rp_i = a_i*(1 - e_i); % Pericentro dell'orbita iniziale

% Parametri orbitali dell' orbita finale
rp_f = a_f*(1 - e_f); % Pericentro dell'orbitale finale

% Prametri orbitali delle orbite di trasferimento
rp_2 = rp_f;  % Pericentro della seconda orbita di trasferimento
rp_1 = rp_i; % Pericentro della prima orbita di trasferimento

ra_1 = ra_t; % Apocentro della prima orbita di trasferimento
ra_2 = ra_1; % Apocentro della seconda orbita di trasferimento

a_1 = (ra_1 + rp_1)/2; % semi-asse maggiore della prima orbita di trasferimento
a_2 = (ra_2 + rp_2)/2; % semi-asse maggiore della seconda orbita di trasferimento

e_1 = (ra_1 - rp_1)/(ra_1 + rp_1); % eccentricità della prima orbita di trasferimento
e_2 = (ra_2 - rp_2)/(ra_2 + rp_2); % eccentricità della seconda orbita di trasferimento

% Calcolo dei costi di trasferimento (delta_v)
DeltaV1 = abs(sqrt(mu*(2/rp_1 - 1/a_1)) - sqrt(mu*(2/rp_i - 1/a_i)));
DeltaV2 = abs(sqrt(mu*(2/ra_2 - 1/a_2)) - sqrt(mu*(2/ra_1 - 1/a_1)));
DeltaV3 = abs(sqrt(mu*(2/rp_f - 1/a_f)) - sqrt(mu*(2/rp_2 - 1/a_2)));

% Calcolo dei tempi di trasferimento (delta_T)
Deltat1 = pi*sqrt(a_1^3/mu);
Deltat2 = pi*sqrt(a_2^3/mu);
DeltaT = Deltat1 + Deltat2;

end

