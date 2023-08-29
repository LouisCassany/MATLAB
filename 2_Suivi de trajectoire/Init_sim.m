clc
clearvars
close all

[r, params] = init_model();
PLOT_VIDEO = false;

controlSaturation = 0.5;
constantControlInput = 0;

%% Init model
theta_0 = 20*pi/180;
phi_0 = 0*pi/180;
dtheta_0 = 10*pi/180;
dphi_0 = 0*pi/180;
psi_0 = 0*pi/180;
x0 = [theta_0; phi_0; dtheta_0; dphi_0; psi_0];

%% Init trajectory
% Génération de trajectoire

% Paramétrage de la trajectoire
s = 21; % nombre de points
a = 40; % Largeur de la trajectoire

% rotation de la trajectoire autour de l'axe des lignes
Omega_0 = 0 * pi / 180;
 
% Centre de la trajectoire, coordonnées polaires
theta_T0 = 60 * pi / 180; 
phi_T0 = 0 * pi / 180;

paramsTrajectoire = [s, a, Omega_0, theta_T0, phi_T0];

Trajectoire = get_trajectory(theta_T0, phi_T0, s, a, Omega_0);

%% Simulation
V_WR = 30/3.6;      % m/s
t_sim = 10;        % s

gain = 2;
epsilon = 0;
d0 = 10;

out = sim('pilote_automatique.slx');

x = out.x.Data;
u = out.u.Data;
time = out.tout(:);
theta = x(:,1);
phi = x(:,2);
dot_theta = x(:,3);
dot_phi = x(:,4);
psi = x(:,5);
error = out.error.Data;
p_kite = r*[(cos(phi).*cos(theta)).'; (sin(phi).*cos(theta)).';-sin(theta).'].';
f=figure;
subplot(1,3,[1 2])
plot_flight_window(r)
hold on
for index = 1:length(Trajectoire)
    plot3(Trajectoire(index).position(1), Trajectoire(index).position(2), Trajectoire(index).position(3), 'b.')
    text(Trajectoire(index).position(1)+2, Trajectoire(index).position(2)+2, Trajectoire(index).position(3)+2,num2str(index))
end
plot3(p_kite(:,1),p_kite(:,2),p_kite(:,3),'m','LineWidth',1,'DisplayName',"Kite's trajectory");
view([-90 30])
rotate3d()
subplot(1,3,3)
hold on
plot(time,theta*180/pi,'DisplayName','\theta (°)')
plot(time,phi*180/pi,'DisplayName','\phi (°)')
plot(time,psi*180/pi,'DisplayName','\psi (°)')
grid on
legend
f.Position = [189 315 1376 420];

% Speed computation

dxdt=diff(p_kite(:,1))./diff(time);
dydt=diff(p_kite(:,2))./diff(time);
dzdt=diff(p_kite(:,3))./diff(time);

v = sqrt(dxdt.^2 + dydt.^2 + dydt.^2);
v=[0; v];

figure
grid on
plot(out.tracking_index)

figure(f)
subplot(1,3,[1 2])
colormap("jet")
a=colorbar;
ylabel(a,'Speed (m/s)','FontSize',12)
clim([min(v), max(v)]);
c = linspace(min(v), max(v), numel(p_kite(:,1)));
surface([p_kite(:,1).';p_kite(:,1).'],[p_kite(:,2).';p_kite(:,2).'],[p_kite(:,3).';p_kite(:,3).'],[v.';v.'],'facecol','no','edgecol','interp','linew',2);
xlabel("X axis")
ylabel("Y axis")
zlabel("Z axis")

%%
if(PLOT_VIDEO)
    name = 'results/essai_lacet_4';
    n=1000;
    frames = animate_kite(n, theta, phi, psi, r);
    speed = 1;
    figure
    plot_flight_window(r)
    movie(frames,1,n/t_sim*speed)
    save(name)
    
    %% Video
    video = VideoWriter(strcat(name,'avi')); %create the video object
    open(video); %open the file for writing
    for ii=1:length(frames) %where N is the number of images
    %   I = imread('the ith image.jpg'); %read the next image
      I = frames(ii);
      writeVideo(video,I); %write the image to file
    end
    close(video); %close the file
end







