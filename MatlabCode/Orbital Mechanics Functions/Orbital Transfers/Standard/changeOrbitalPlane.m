function [deltaV, om_f, theta, alpha] = ...
    changeOrbitalPlane(a, e, i_i, OM_i, om_i, i_f, OM_f, mu)

% Cambio di piano dell'orbita
% 
% [DeltaV, om_f, theta] = ...
%   changeOrbitalPlane(a, e, i_i, OM_i, om_i, i_f, OM_f, mu)
% 
% ---------------------------------------------------------------------------
% Argomenti di Input:
% a                 [1x1]   semi-asse maggiore                   [km]
% e                 [1x1]   eccentricità                         [-]
% i_i               [1x1]   inclinazione iniziale                [rad]
% OM_i              [1x1]   RAAN iniziale                        [rad]
% om_i              [1x1]   anomalia del pericentro iniziale     [rad]
% i_f               [1x1]   semi-asse maggiore                   [rad]
% OM_f              [1x1]   RAAN finale                          [rad]
% mu                [1x1]   costante gravitazionale              [km^3/s^2]
% 
% ---------------------------------------------------------------------------
% Argomenti di Output:
% DeltaV            [1x1]   costo della manovra                  [km]
% om_f              [1x1]   anomalia del pericentro finale       [rad]
% theta             [1x1]   anomalia vera della manovra          [rad]
% alpha             [1x1]   anomalia formata dai due piano       [rad]
%
% ---------------------------------------------------------------------------

% Risoluzione del triangolo sferico

deltaOM = OM_f - OM_i;
deltai = i_f - i_i;

if((deltaOM) > 0)
    if(deltai > 0)
        alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(deltaOM));
        sin_u_i = sin(i_f)*(sin(deltaOM)/sin(alpha));
        cos_u_i = (-cos(i_f) + cos(alpha)*cos(i_i))/(sin(alpha)*sin(i_i));
        u_i = atan2(sin_u_i, cos_u_i);
        sin_u_f = sin(i_i)*(sin(deltaOM)/sin(alpha));
        cos_u_f = (cos(i_i) - cos(alpha)*cos(i_f))/(sin(alpha)*sin(i_f));
        u_f = atan2(sin_u_f, cos_u_f);
        theta1 = u_i - om_i;
        om_f = u_f - theta1;
    else
        alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(deltaOM));
        sin_u_i = sin(i_f)*(sin(deltaOM)/sin(alpha));
        cos_u_i = (cos(i_f) - cos(alpha)*cos(i_i))/(sin(alpha)*sin(i_i));
        u_i = atan2(sin_u_i, cos_u_i);
        sin_u_f = sin(i_i)*(sin(deltaOM)/sin(alpha));
        cos_u_f = (-cos(i_i) + cos(alpha)*cos(i_f))/(sin(alpha)*sin(i_f));
        u_f = atan2(sin_u_f, cos_u_f);
        theta1 = 2*pi - u_i - om_i;
        om_f = 2*pi - u_f - theta1;
    end
else
    deltaOM = abs(deltaOM);
    if(deltai > 0)
        alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(deltaOM));
        sin_u_i = sin(i_f)*(sin(deltaOM)/sin(alpha));
        cos_u_i = (-cos(i_f) + cos(alpha)*cos(i_i))/(sin(alpha)*sin(i_i));
        u_i = atan2(sin_u_i, cos_u_i);
        sin_u_f = sin(i_i)*(sin(deltaOM)/sin(alpha));
        cos_u_f = (cos(i_i) - cos(alpha)*cos(i_f))/(sin(alpha)*sin(i_f));
        u_f = atan2(sin_u_f, cos_u_f);
        theta1 = 2*pi - u_i - om_i;
        om_f = 2*pi - u_f - theta1;
    else
        alpha = acos(cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(deltaOM));
        sin_u_i = sin(i_f)*(sin(deltaOM)/sin(alpha));
        cos_u_i = (cos(i_f) - cos(alpha)*cos(i_i))/(sin(alpha)*sin(i_i));
        u_i = atan2(sin_u_i, cos_u_i);
        sin_u_f = sin(i_i)*(sin(deltaOM)/sin(alpha));
        cos_u_f = (-cos(i_i) + cos(alpha)*cos(i_f))/(sin(alpha)*sin(i_f));
        u_f = atan2(sin_u_f, cos_u_f);
        theta1 = u_i - om_i;
        om_f = u_f - theta1;
    end
end

% Calcolo del costo della manovra
p = a*(1 - e^2); % semilato retto
v_theta1 = sqrt(mu/p)*(1 + e*cos(theta1)); % velocità nel primo punto di cambio piano
theta2 = theta1 + pi;
v_theta2 = sqrt(mu/p)*(1 + e*cos(theta2)); % velocità nel secondo punto di cambio piano

% Costo del cambio di piano
if(v_theta1 > v_theta2)
    deltaV = 2*v_theta2*sin(alpha/2); 
    theta = theta2;
else
    deltaV = 2*v_theta1*sin(alpha/2);
    theta = theta1;
end

end

