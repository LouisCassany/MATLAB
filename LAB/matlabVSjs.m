clc
clear all
close all

r = 50;
theta = 45*pi/180;
phi = 0*pi/180;
n = 20;
a = 20;
psi = 0*pi/180;

Trajectoire = get_trajectory(theta, phi, n, a, psi);

figure;
plot_flight_window(r)
hold on
grid on
axis equal

for P = Trajectoire
    plot3(P.position(1), P.position(2), P.position(3),'.b')
end

rotate3d()