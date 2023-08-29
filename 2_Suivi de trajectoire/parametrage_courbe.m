clc
clear all
close all

r = 50;

Mx = @(alpha) [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
My = @(alpha) [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
Mz = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

P = [[0 0] ; [-20 20] ; [0 0]];

f=figure;
rotate3d(f)
f.Position = [361 66 803 658];
plot_flight_window(r)
hold on
plot3(P(1,:), P(2,:), P(3,:),'.b',MarkerSize=50)
view([270 30])

theta_T0 = 70*pi/180;
phi_T0 = 0*pi/180;
psi_T0 = 0*pi/180;

Pc = [r*cos(theta_T0)*cos(phi_T0); r*cos(theta_T0)*sin(phi_T0);-r*sin(theta_T0)];

% Matrice de passage de R_WR -> R_k0 + rotation
WR_M_k0 = Mz(phi_T0)*My(theta_T0-pi/2)*Mz(psi_T0);

P2 =  WR_M_k0*P + Pc;

plot3(P2(1,:), P2(2,:), P2(3,:),'.r', MarkerSize=50)

P3 = P2./vecnorm(P2)*50;

plot3(P3(1,:), P3(2,:), P3(3,:),'.g', MarkerSize=50)

v1 = P3(:,2);
v2 = P3(:,1);

plot_vector([0;0;0], v1,'b','v1')
plot_vector([0;0;0], v2,'r','v2')

signed_angle = signed_angle_between_vectors(v1, v2, [0;0;-10])*180/pi











