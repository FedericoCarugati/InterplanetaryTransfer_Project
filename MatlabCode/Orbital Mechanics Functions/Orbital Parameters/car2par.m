function [a, e, i, OM, om, theta] = car2par(r,v,mu)

% Trasformazione dalle coordinate cartesiane ai parametri Kepleriani
%
% [a, e, i, OM, om, th] = car2par(rr, vv, mu)
%
% --------------------------------------------------------------------------
% Argomenti di Input:
% rr               [3x1]   vettore posizione                     [km]
% vv               [3x1]   vettore velocità                      [km/s]
% mu               [3x1]   costante gravitazionale               [km^3/s^2]
%
% --------------------------------------------------------------------------
% Argomenti di Output:
% a                [1x1]   semi-asse maggiore                    [km]
% e                [1x1]   eccentricità                          [-]
% i                [1x1]   inclinazione                          [rad]
% OM               [1x1]   ascensione retta del nodo ascendente  [rad]
% om               [1x1]   anomalia del pericentro               [rad]
% th               [1x1]   anomalia vera                         [rad]
%
% --------------------------------------------------------------------------

r_norm = norm(r);
v_norm = norm(v);

k = [0 0 1]';

% E = 0.5*(v_norm)^2 - mu/r_norm; Enenrgia specifica dell'orbita
a = (2/r_norm - v_norm^2/mu)^(-1); % Semiasse maggiore dell'orbita

% Vettore momento angolare
h = cross(r,v);
h_norm = norm(h);

% Vettore eccentricità
e_vect = (cross(v,h))/mu - r/r_norm;
e = norm(e_vect);

% Inclinazione dell'orbita
i = acos(h(3)/h_norm);

% Nodo ascendent - Linea dei nodi
N = (cross(k,h)/norm(cross(k,h)));

% Ascensione retta del nodo ascendente
if(N(2) >= 0)
    OM = acos(N(1));
else
    OM = 2*pi - acos(N(1));
end

% Anomalia del pericentro
if(e_vect(3) >= 0)
    om = acos(dot(N,e_vect)/e);
else
    om = 2*pi - acos(dot(N,e_vect)/e);
end

v_r = (dot(v,r))/r_norm;

% Anomalia vera
if(v_r >= 0)
    theta = acos(dot(r,e_vect)/(r_norm*e));
else
    theta = 2*pi - acos(dot(r,e_vect)/(r_norm*e));
end

end

