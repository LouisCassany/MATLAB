function [frames] = animate_kite(n,theta,phi,psi,r)
% ANIMATE_KITE(n,theta,phi,psi) creates a n frames animation of the
% kite's trajectory. r is the tether length and theta, phi, psi are column 
% vectors of the elevation, azimuth, yaw of the kite in radians.
%
% You may use 'movie(frames)' to play the animation.

p_kite = r*[cos(phi).*cos(theta) sin(phi).*cos(theta) -sin(theta)];

WR_M_k0 = @(phi_arg, theta_arg) [cos(phi_arg)*sin(theta_arg) -sin(phi_arg) -cos(phi_arg)*cos(theta_arg); sin(phi_arg)*sin(theta_arg) cos(phi_arg) -sin(phi_arg)*cos(theta_arg); cos(theta_arg) 0 sin(theta_arg)];
k0_M_b = @(psi_arg) [-cos(psi_arg) sin(psi_arg) 0; -sin(psi_arg) -cos(psi_arg) 0; 0 0 1];

ind = round(linspace(1,length(p_kite(:,1)),n));

f = figure;
plot_flight_window(r)
hold on

view([270 30])
axis tight manual

frames(n) = struct('cdata',[],'colormap',[]);

h = waitbar(0,'Animation in progress');
f.Visible = 'off';

for j = 1:n
    t=ind(j);
    p1 = WR_M_k0(phi(t),theta(t))*k0_M_b(psi(t))*[4;0;0] + p_kite(t,:).';
    p2 = WR_M_k0(phi(t),theta(t))*k0_M_b(psi(t))*[-2;2;0] + p_kite(t,:).';
    p3 = WR_M_k0(phi(t),theta(t))*k0_M_b(psi(t))*[-2;-2;0] + p_kite(t,:).';
    % p1 = WR_M_k0(phi(t),theta(t))*k0_M_b(psi(t))*[2;4;0] + p_kite(t,:).';
    % p2 = WR_M_k0(phi(t),theta(t))*k0_M_b(psi(t))*[2;-4;0] + p_kite(t,:).';
    % p3 = WR_M_k0(phi(t),theta(t))*k0_M_b(psi(t))*[-2;-4;0] + p_kite(t,:).';
    % p4 = WR_M_k0(phi(t),theta(t))*k0_M_b(psi(t))*[-2;4;0] + p_kite(t,:).';
    p = [p1 p2 p3];
    b=fill3(p(1,:),p(2,:),p(3,:),'g');
    view([270 30])
    drawnow
    frames(j) = getframe;
    waitbar(j/n);
    b.Visible = 'off';
end
close(f)
close(h)
end

