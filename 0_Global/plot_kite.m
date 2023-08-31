function plot_kite(theta,phi,psi,r, size)
% PLOT_KITE(theta,phi,psi,r)
%
% Plot the kite at the given position.
%
% theta : elevation (rad)
% phi   : phi (rad)  
% psi   : psi (rad)
% r     : tether length (m)
% size  : kite size (m)

p_kite = r*[cos(phi).*cos(theta) sin(phi).*cos(theta) -sin(theta)];

WR_M_k0 = @(phi_arg, theta_arg) [cos(phi_arg)*sin(theta_arg) -sin(phi_arg) -cos(phi_arg)*cos(theta_arg); sin(phi_arg)*sin(theta_arg) cos(phi_arg) -sin(phi_arg)*cos(theta_arg); cos(theta_arg) 0 sin(theta_arg)];
k0_M_b = @(psi_arg) [-cos(psi_arg) sin(psi_arg) 0; -sin(psi_arg) -cos(psi_arg) 0; 0 0 1];

p1 = WR_M_k0(phi,theta)*k0_M_b(psi)*[size/2;-size/2;0] + p_kite.';
p2 = WR_M_k0(phi,theta)*k0_M_b(psi)*[size/2;size/2;0] + p_kite.';
p3 = WR_M_k0(phi,theta)*k0_M_b(psi)*[-size/2;size;0] + p_kite.';
p4 = WR_M_k0(phi,theta)*k0_M_b(psi)*[-size/2;-size;0] + p_kite.';

p = [p1 p2 p3 p4];
% b=fill3(p(1,:),p(2,:),p(3,:),'g');

plot3(p_kite(1), p_kite(2), p_kite(3), 'r.', MarkerSize=20)

plot_vector(p_kite, WR_M_k0(phi, theta)*[10; 0; 0],'k','x_{k0}')
plot_vector(p_kite, WR_M_k0(phi, theta)*[0; 10; 0],'k','y_{k0}')
plot_vector(p_kite, WR_M_k0(phi, theta)*[0; 0; 10],'k','z_{k0}')

plot_vector(p_kite, WR_M_k0(phi, theta)*k0_M_b(psi)*[20; 0; 0],'b','x_{b}')
plot_vector(p_kite, WR_M_k0(phi, theta)*k0_M_b(psi)*[0; 20; 0],'b','y_{b}')
plot_vector(p_kite, WR_M_k0(phi, theta)*k0_M_b(psi)*[0; 0; 20],'b','z_{b}')


rotate3d










