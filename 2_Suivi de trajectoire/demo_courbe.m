clc
clear all
close all

[r,~] = init_model();

% Centre de la trajectoire, coordonnées polaires
theta_T0 = 50*pi/180;
phi_T0 = 30*pi/180;

% Centre de la trajectoire, coordonnées cartésiennes (dans R_WR)
Ch = [r*cos(theta_T0)*cos(phi_T0); r*cos(theta_T0)*sin(phi_T0);-r*sin(theta_T0)];

% Matrice de passages
Mx = @(alpha) [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
My = @(alpha) [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
Mz = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

% Matrice de passage de R_WR -> R_k0
WR_M_k0 = Mz(phi_T0)*My(theta_T0-pi/2);

% Paramétrage de la trajectoire
s = linspace(0,2*pi,100);       % Discretisation
a = 20;                         % Largeur de la trajectoire
Omega_0 = 0*pi/180;             % rotation de la trajectoire autour de l'axe des lignes

%% 1. Coordonnées des points de la trajectoire dans k0;
yTdes = a*(sin(s))./(1+cos(s).^2);
xTdes = a*(sin(s).*cos(s))./(1+cos(s).^2);
k0_Tdes = [xTdes;yTdes;zeros(1,length(xTdes))];

% plot
f=figure;
rotate3d(f)
f.Position = [361 66 803 658];
subplot(221)
plot_flight_window(r)
hold on
plot3(k0_Tdes(1,:),k0_Tdes(2,:),k0_Tdes(3,:))
plot3(Ch(1),Ch(2),Ch(3),'g.','MarkerSize',20');
title('1. Trajectoire dans R_{k0}')
view([270 30])

%% 2. Passage dans R_WR
WR_Tdes = [];
for P=k0_Tdes
    WR_Tdes = [WR_Tdes WR_M_k0*P+Ch];
end

% plot
subplot(222)
plot_flight_window(r)
hold on
plot3(WR_Tdes(1,:),WR_Tdes(2,:),WR_Tdes(3,:))
plot3(Ch(1),Ch(2),Ch(3),'g.','MarkerSize',20');
title('2. Projection dans R_{WR}')
view([270 30])


%% 3. Rotation de la trajectoire autour de l'axe des lignes

% Matrice de passage de R_WR -> R_k0 + rotation
WR_M_k0 = Mz(phi_T0)*My(theta_T0-pi/2)*Mz(Omega_0);

yTdes = a*(sin(s))./(1+cos(s).^2);
xTdes = a*(sin(s).*cos(s))./(1+cos(s).^2);
k0_Tdes = [xTdes;yTdes;zeros(1,length(xTdes))];

WR_Tdes = [];
for P=k0_Tdes
    WR_Tdes = [WR_Tdes WR_M_k0*P + Ch];
end

% plot
subplot(223)
plot_flight_window(r)
hold on
plot3(WR_Tdes(1,:),WR_Tdes(2,:),WR_Tdes(3,:))
plot3(Ch(1),Ch(2),Ch(3),'g.','MarkerSize',20');
title("3. Rotation d'Omega_0")
view([270 30])

%% 4. Projection sur la fen�tre de vol

WR_Tdes_norm = [];
for P=WR_Tdes
    WR_Tdes_norm = [WR_Tdes_norm P./norm(P)*r];
end

subplot(224)
plot_flight_window(r)
hold on
plot3(WR_Tdes_norm(1,:),WR_Tdes_norm(2,:),WR_Tdes_norm(3,:),'b')
plot3(Ch(1),Ch(2),Ch(3),'g.','MarkerSize',20');
title('4. Projection sur la fenêtre de vol')
view([270 30])




