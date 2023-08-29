clc
clear all
close all

theta_T0 = 30*pi/180;
phi_T0 = 0*pi/180;
n = 19;
a = 50;
Omega_T0 = 0;
d0 = 15;
K_p = 1;
r = 50;
last_target_index = 1;
theta = 18*pi/180;
phi = 12*pi/180;
psi = 0;

p_kite = r*[(cos(phi).*cos(theta)).'; (sin(phi).*cos(theta)).';-sin(theta).'].';


[control_signal, alpha, alpha_des, new_target_index, error, Trajectory] = pilote_auto(...
    theta_T0, phi_T0, n, a, Omega_T0, ...
    d0, K_p, r,...
    1, ...
    theta, phi, psi);

target = new_target_index

f = figure;
plot_flight_window(r)
hold on
plot3(p_kite(1), p_kite(2), p_kite(3), 'm.', "MarkerSize", 30)
grid on
for P = Trajectory
    plot3(P.position(1), P.position(2), P.position(3),'r.', "MarkerSize", 20)
    text(P.position(1), P.position(2), P.position(3),num2str(P.index))
end
view([90 0])
rotate3d

f.Position = [288, 191, 1146, 726];

% for i = 1:10
%     target
%     [control_signal, alpha, alpha_des, new_target_index, error, Trajectory] = pilote_auto(...
%     theta_T0, phi_T0, n, a, Omega_T0, ...
%     d0, K_p, r,...
%     target, ...
%     theta, phi, psi);
% 
%     target = new_target_index
% end

plot_vector(p_kite, Trajectory(new_target_index).position - p_kite.', 'b')




