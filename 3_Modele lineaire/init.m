clc
clear all
close all

%% Load and plot the data
load("data_pilote_auto.mat");

figure
subplot(211)
hold on
grid on
plot(data.time, data.theta,"DisplayName","Élévation")
plot(data.time, data.phi,"DisplayName","Azimuth")
plot(data.time, data.psi,"DisplayName","Lacet")
xlabel('Time (s)');
ylabel('Angle (°)')
legend
title("États kite")
subplot(212)
plot(data.time, data.delta)
grid on
xlabel('Time (s)');
ylabel('Différentiel (mm)')
title("Signal de commande")

%% FFT

Fs = 5; % Sampling frequency (Hz)
t = 0:1/Fs:data.time(end);
nfft = 1024;

u = timeseries(data.delta, data.time);
u = resample(u, t);
y = timeseries(data.psi, data.time);
y = resample(y, t);

U = fft(u.Data,nfft);
U = U(1:nfft/2);
mu = abs(U);

Y = fft(y.Data,nfft);
Y = Y(1:nfft/2);
my = abs(Y);
f = (0:nfft/2-1)*Fs/nfft;

figure
subplot(2,4,[1 2 5 6])
hold on
grid on
plot(u)
plot(y)
legend("Commande","Lacet")
xlabel("Time (s)")

subplot(2,4,[3 4 7 8])
hold on
plot(f,mu);
plot(f,my);
xlabel("Frequency (Hz)")
ylabel("Power")
grid on
legend("Commande","Lacet")

%% Linearization

load("modelABCD.mat");

% % Constants parameters
% [~,params]= init_model();
% r = params(1);
% m = params(2);
% g_k = params(3);
% M_k = params(4);
% A_k = params(5);
% rho_air = params(6);
% g = params(7);
% K_i = params(8);
% alpha_i0 = params(9);
% C_D = params(10);
% C_L = params(11);
% V_WR = 30/3.6;          % m/s
% 
% % Transformation matrices
% Mx = @(alpha) [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
% My = @(alpha) [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
% Mz = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
% 
% syms theta phi dtheta dphi psi delta epsilon real
% 
% WR_M_k0 = Mz(phi)*My(theta-pi/2); % from Rk0 to RWR
% V_WR = [V_WR;0;0];
% V_k_k0 = [-r*dtheta;r*dphi*cos(theta);0]; % Kite's speed in Rk0
% V_a_WR = V_WR-WR_M_k0*V_k_k0;
% modV_a = norm(V_a_WR);
% 
% k0_M_b = Mz(pi+psi);
% xb_k0 = k0_M_b*[1;0;0];
% yb_k0 = k0_M_b*[0;1;0];
% zb_k0 = k0_M_b*[0;0;1];
% xa_WR = -V_a_WR/modV_a;
% xa_k0 = WR_M_k0.'*xa_WR;
% 
% alpha_0 = pi/2 - acos(dot(zb_k0,xa_k0)/(norm(zb_k0)*norm(xa_k0)));
% alpha_i = K_i*epsilon+alpha_i0;
% alpha = alpha_0 + alpha_i;
% 
% lift = 0.5*rho_air*C_L*alpha*A_k*modV_a*modV_a;
% drag = 0.5*rho_air*C_D*alpha*A_k*modV_a*modV_a;
% 
% x_lift_k0=cross(-xa_k0,yb_k0)/norm(cross(-xa_k0,yb_k0));
% x_drag_k0= -xa_k0;
% 
% F_a_k0 = lift*(x_lift_k0)+drag*(x_drag_k0);
% 
% P_k0 = [m*g*cos(theta); 0; m*g*sin(theta)];
% 
% F_ext_k0 = P_k0 + F_a_k0;
% 
% V_a_Rb=k0_M_b'*WR_M_k0'*V_a_WR;
% va = -V_a_Rb(1);
% 
% ddtheta = (-F_ext_k0(1)/(r*m))-sin(theta)*cos(theta)*dphi^2;
% ddphi = (F_ext_k0(2)/(r*m*cos(theta)))+ 2*tan(theta)*dtheta*dphi;
% dpsi = g_k*va*delta+(M_k*((cos(theta)*sin(psi))/va))-dphi*sin(theta);
% 
% dot_x = [dtheta; dphi; ddtheta; ddphi; dpsi].';
% 
% dot_x = eval(dot_x);
% 
% 
% x = [theta phi dtheta dphi psi].';
% u = delta;
% y = [theta; phi; psi];
% 
% A = jacobian(dot_x,x);
% B = jacobian(dot_x,u);
% C = jacobian(y,x);
% D = jacobian(y,u);
% 
% save('modelABCD','A','B','C','D');

theta = 60*pi/180;
phi = 0;
psi = 10*pi/180;
dtheta = 0;
dphi = 0;
delta = 0;
epsilon = 0;

Aeval = eval(A);
Beval = eval(B);
Ceval = eval(C);
Deval = eval(D);

sys = ss(Aeval, Beval, Ceval, Deval);

simin = [data.time, data.delta];




