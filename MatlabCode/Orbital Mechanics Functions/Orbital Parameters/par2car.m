function [r, v] = par2car(a, e, i, omega, w, theta, mu)

% Trasformazione dai parametri Kepleriani alle coordinate cartesiane
%
% -----------------------------------------------------------------------------
% Imput:
% a               [1x1]  semi asse maggiore                          [km]
% e               [1x1]  eccentricità                                [-]
% i               [1x1]  inclinazione                                [rad]
% omega           [1x1]  acensione retta del nodo ascendente         [rad]
% w               [1x1]  anomalia del pericentro                     [rad]
% theta           [1x1]  anomalia vera                               [rad]
% mu              [1x1]  parametro gravitazionale                    [km^3/s^2]
%
% ------------------------------------------------------------------------------
%Outputs:
% r               [3x1]  vettore posizione                           [km]
% v               [3x1]  vettore velocità                            [km/s]
%
% ------------------------------------------------------------------------------

p = a*(1 - e^2);  % Semilato retto
r_abs = p/(1 + e*cos(theta)); % Modulo vettore posizione

% Matrici di rotazione:
R_om = [cos(omega)  sin(omega)      0;
        -sin(omega) cos(omega)      0;
            0           0           1];

R_i = [     1           0           0;
            0         cos(i)      sin(i);
            0        -sin(i)      cos(i)];

R_w = [  cos(w)       sin(w)        0;
        -sin(w)       cos(w)        0;
           0            0           1];

T = R_w*R_i*R_om;

% Vettori di stato nel sistema perifocale
r_pf = r_abs*[cos(theta); sin(theta); 0];
v_pf = sqrt(mu/p)*[-sin(theta); e+cos(theta); 0];

% Vettori di stato nel sistema di rif. celeste
r = T'*r_pf;
v = T'*v_pf;

end


