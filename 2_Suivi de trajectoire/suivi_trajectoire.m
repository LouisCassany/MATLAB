clc
clearvars
close all

[r, ~] = init_model();

WR_M_k0 = @(phi_arg, theta_arg) [cos(phi_arg)*sin(theta_arg) -sin(phi_arg) -cos(phi_arg)*cos(theta_arg); sin(phi_arg)*sin(theta_arg) cos(phi_arg) -sin(phi_arg)*cos(theta_arg); cos(theta_arg) 0 sin(theta_arg)];
k0_M_b = @(psi_arg) [-cos(psi_arg) sin(psi_arg) 0; -sin(psi_arg) -cos(psi_arg) 0; 0 0 1];


% Position du kite pour test
theta = 60 * pi / 180;
phi = -40 * pi / 180;
psi = 50 * pi / 180;
P_kite = r * [(cos(phi) .* cos(theta)).'; (sin(phi) .* cos(theta)).'; -sin(theta).'];

V_kite(1) = -6;
V_kite(2) = 2;
V_kite(3) = 0;

% Génération de trajectoire

% Paramétrage de la trajectoire
s = 20; % nombre de points
a = 30; % Largeur de la trajectoire

% rotation de la trajectoire autour de l'axe des lignes
Omega_0 = 0 * pi / 180;
 
% Centre de la trajectoire, coordonnées polaires
theta_T0 = 60 * pi / 180; 
phi_T0 = 0 * pi / 180;

Trajectoire = get_trajectory(theta_T0, phi_T0, s, a, Omega_0);

f= figure;
hold on
rotate3d()
view(-90, 90);
grid on
ax = gca;
ax.XDir = 'Reverse';
ax.ZDir = 'Reverse';
axis equal

% Plot trajectoire
for index = 1:length(Trajectoire)
    plot3(Trajectoire(index).position(1), Trajectoire(index).position(2), Trajectoire(index).position(3), 'b.')
    text(Trajectoire(index).position(1), Trajectoire(index).position(2), Trajectoire(index).position(3), num2str(index))
end


%% Plots kite et P_i
plot_flight_window(r)
hold on

% Asservissement en lacet
d0 = 2;
last_index = 1;
[alpha, alpha_des, new_index] = asservissement_lacet(P_kite, Trajectoire, d0, last_index, theta, phi, psi);








