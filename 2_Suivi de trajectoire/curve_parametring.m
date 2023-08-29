
clc
clear all
close all

% Matrice de rotations
Mx = @(alpha) [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
My = @(alpha) [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
Mz = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

r = 50;

% Centre de la trajectoire, coordonnées polaires
theta_T0 = 30*pi/180;
phi_T0 = 20*pi/180;

% Paramétrage de la trajectoire
s = linspace(0,2*pi,20);       % Discretisation
a = 20;                         % Largeur de la trajectoire
Omega_0 =30*pi/180;             % rotation de la trajectoire autour de l'axe des lignes

% Centre de la trajectoire, coordonnées cartésiennes (dans R_WR)
Ch = [r*cos(theta_T0)*cos(phi_T0); r*cos(theta_T0)*sin(phi_T0);-r*sin(theta_T0)];

yTdes = a*(sin(s))./(1+cos(s).^2);
xTdes = a*(sin(s).*cos(s))./(1+cos(s).^2);
zTdes = zeros(1,length(xTdes));
k0_Tdes = [xTdes;yTdes;zTdes];


ka_Tdes = [0, a, 0].';

f=figure
% subplot(121)
plot_flight_window(r)
hold on
grid on
plot(k0_Tdes(1,:), k0_Tdes(2,:),'.k')
plot(ka_Tdes(1), ka_Tdes(2),'.b','MarkerSize',30)
plot_vector([-15,-20,0],[10,0,0],'r','X')
plot_vector([-15,-20,0],[0,10,0],'g','Y')
view([90 0])
title('R_{WR}', 'Interpreter','tex')
rotate3d


f.Position = [286, 300, 1183, 400];

WR_M_k0 = Mz(phi_T0)*My(theta_T0-pi/2)*Mz(Omega_0);

WR_Tdes = (WR_M_k0*k0_Tdes+Ch);
WR_Tdes = WR_Tdes./vecnorm(WR_Tdes)*r;

ka_Tdes = (WR_M_k0*ka_Tdes+Ch);
ka_Tdes = ka_Tdes./vecnorm(ka_Tdes)*r;





% subplot(122)
plot_flight_window(r)
hold on
grid on
plot3(WR_Tdes(1,:),WR_Tdes(2,:),WR_Tdes(3,:),'.k')
plot3(ka_Tdes(1,:),ka_Tdes(2,:),ka_Tdes(3,:),'.b','MarkerSize',30)
view([90 0])





