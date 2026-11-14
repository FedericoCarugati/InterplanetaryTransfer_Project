function plotorbit_hyperbolic1(a, e, i, OM, om, delta_th, dim, mu)

% 3D Plot dell'orbita - linea continua
%
% plotorbit_hyperbolic1(a, e, i, OM, om, delta_th, dim, mu)
%
% -------------------------------------------------------------------------------
%
% Argomenti di Input
% a                    [1x1]    semi-asse maggiore                     [km]
% e                    [1x1]    eccentricità                           [-]
% i                    [1x1]    inclinazione                           [rad]
% OM                   [1x1]    ascensione retta del nodo ascendente   [rad]
% om                   [1x1]    anomalia del pericentro                [rad]
% delta_th             [1x1]    step size dell'anomalia vera           [rad]
% dim                  [1x1]    dimensione caratteristica              [km]
% mu                   [1x1]    costante gravitazionale                [km^3/s^2]
% -------------------------------------------------------------------------------

p = a*(1 - e^2); % semilato retto dell'orbita

th0 = 0;
thf = acos((p/dim - 1)/e);
n = round((thf - th0)/delta_th);
theta = linspace(th0, thf, n);

% Inizializzazione dei vettori cartesiani
r = zeros(3,n);
v = zeros(3,n);

% Costruzione dei vettori cartesiani di velocità e posizione
for(j = 1:n)
    [r_i, v_i] = par2car(a, e , i, OM, om, theta(j), mu);
    r(:,j) = r_i;
    v(:,j) = v_i;
end

% Plot dell'orbita
plot3(r(1,:),r(2,:),r(3,:),'LineWidth',1.0)
hold on
grid on
axis equal
xlabel('X [km]')
ylabel('Y [km]')
zlabel('Z [km]')
% scatter3(r(1,1), r(2,1), r(3,1), 30, 'filled')
% scatter3(r(1,end), r(2,end), r(3,end), 30, 'filled')

end

