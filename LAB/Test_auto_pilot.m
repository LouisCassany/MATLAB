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
kite_size = 3;

theta = 20*pi/180;
phi = 20*pi/180;
psi = 40*pi/180;

p_kite = r*[(cos(phi).*cos(theta)).'; (sin(phi).*cos(theta)).';-sin(theta).'].';

WR_M_k0 = @(phi_arg, theta_arg) [cos(phi_arg)*sin(theta_arg) -sin(phi_arg) -cos(phi_arg)*cos(theta_arg); sin(phi_arg)*sin(theta_arg) cos(phi_arg) -sin(phi_arg)*cos(theta_arg); cos(theta_arg) 0 sin(theta_arg)];
k0_M_b = @(psi_arg) [-cos(psi_arg) sin(psi_arg) 0; -sin(psi_arg) -cos(psi_arg) 0; 0 0 1];

Trajectory = get_trajectory(theta_T0, phi_T0, n, a, Omega_T0);

[control_signal, alpha, alpha_des, target_index, next_target_index, error] = pilote_auto(...
    theta_T0, phi_T0, n, a, Omega_T0, ...
    d0, K_p, r,...
    1, ...
    theta, phi, psi);

f = figure;
plot_flight_window(r)
hold on
plot_kite(theta, phi, psi, r, kite_size)
grid on
for P = Trajectory
    plot3(P.position(1), P.position(2), P.position(3),'r.', "MarkerSize", 20)
    text(P.position(1), P.position(2), P.position(3),num2str(P.index))
end

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

plot_vector(p_kite, Trajectory(target_index).position - p_kite.', 'b')
plot_vector(p_kite, Trajectory(next_target_index).position - p_kite.', 'r')

WR_V_next = WR_M_k0(phi, theta)*[2.6827; 9.6585; 1.4591];

plot_vector(p_kite, WR_V_next,'g','V_{next}');

view([90 0])
rotate3d




