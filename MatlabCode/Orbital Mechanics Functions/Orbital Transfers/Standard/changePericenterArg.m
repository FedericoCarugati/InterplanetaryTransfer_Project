function [dv, theta_i, theta_f] = changePericenterArg(a, e, wi, wf, mu)

% Manovra di cambio di argomento del pericentro
%
% --------------------------------------------------------------------------
% Argomenti di input:
% a               [1x1]  semi-asse maggiore                          [km]
% e               [1x1]  eccentricità                                [-]
% wi              [1x1]  anomalia del pericentro iniziale            [rad]
% wf              [1x1]  anomalia del pericentro finale              [rad]
% mu              [1x1]  parametro gravitazionale                    [km^3/s^2]
%
% --------------------------------------------------------------------------
% Argomenti di output:
% dv              [1x1]  impulso di manovra                          [km/s]
% theta_i         [2x1]  anomalia vera iniziale                      [rad]
% theta_f         [2x1]  anomalia vera finale                        [rad]

dw = wf - wi; % delta di anomalie di pericentro
p = a*(1 - e^2);  % semilato retto

% Calcolo delle anomalie vere dove la manovra può essere effettuata
th1i = dw/2;
th2i = pi + dw/2;

th1f = 2*pi - dw/2;
th2f = pi - dw/2;

theta_i = [th1i; th2i];
theta_f = [th1f; th2f];

% Costo della manovra
dv = abs(2*sqrt(mu/p)*e*sin(dw/2));

end


