clc
clearvars
close all

%% General tweaks
fontSize = 13;
show_animation = true;

%% Init model
theta_0 = 65*pi/180;
phi_0 = -10*pi/180;
dtheta_0 = 0*pi/180;
dphi_0 = 0*pi/180;
psi_0 = 0*pi/180;
x0 = [theta_0; phi_0; dtheta_0; dphi_0; psi_0];
[r,params] = init_model();

V_WR = 30/3.6;      % km/h

% stepSize = 0.1; % s
t_sim = 30;     % s

out = sim('kite');

%% plot results
x = out.x.Data;
u = out.u.Data;
time = out.tout(:);
theta = x(:,1);
phi = x(:,2);
dot_theta = x(:,3);
dot_phi = x(:,4);
psi = x(:,5);
p_kite = r*[(cos(phi).*cos(theta)).'; (sin(phi).*cos(theta)).';-sin(theta).'].';

f=figure;
f.Position = [f.Position(1) f.Position(2) f.Position(3)*2 f.Position(4)];
subplot(1,2,1);
plot(time, theta*180/pi,'b');
hold on
plot(time, dot_theta*180/pi,'b--');
plot(time, phi*180/pi,'r');
plot(time, dot_phi*180/pi,'r--');
plot(time, psi*180/pi);
grid on
xlabel('Time (s)','FontSize',fontSize)
ylabel('System states','FontSize',fontSize)
l=legend('$\theta(^{\circ})$','$\dot{\theta}$ $(^{\circ}/s)$','$\phi$ $(^{\circ})$','$\dot{\phi}$ $(^{\circ}/s)$','$\psi$ (${^\circ})$');
l.Location = "best";
l.FontSize = fontSize;
l.Interpreter = "latex";

% truncate trajectory vector if theta < 0
vals = find(theta < 0);
if ~isempty(vals)
    index = vals(1);
    p_kite = p_kite(1:index-1,:);
end

% Plot kite trajectory
subplot(1,2,2)
plot_flight_window(r);
hold on
view([270 30])
plot3(p_kite(:,1),p_kite(:,2),p_kite(:,3),'m--','LineWidth',2,'DisplayName',"Kite's trajectory");
hold on
plot3(p_kite(1,1),p_kite(1,2),p_kite(1,3),'g.','MarkerSize',30,'DisplayName',"Initial position");
plot3(p_kite(end,1),p_kite(end,2),p_kite(end,3),'r.','MarkerSize',30,'DisplayName',"Final position");
grid on
rotate3d(f)
l=legend;
l.Location = "northeast";
l.FontSize = fontSize;
title("Kite trajectory in flight window",'FontSize',fontSize)

%%
if(show_animation)
    n=100;
    frames = animate_kite(n, theta, phi, psi, r);
    speed = 4;
    figure
    movie(frames,1,n/t_sim*speed)
end




