function [alpha, alpha_des, new_index, error] = asservissement_lacet(P_kite, Trajectoire, d0, last_index, theta, phi, psi)

WR_M_k0 = @(phi_arg, theta_arg) [cos(phi_arg)*sin(theta_arg) -sin(phi_arg) -cos(phi_arg)*cos(theta_arg); sin(phi_arg)*sin(theta_arg) cos(phi_arg) -sin(phi_arg)*cos(theta_arg); cos(theta_arg) 0 sin(theta_arg)];
k0_M_b = @(psi_arg) [-cos(psi_arg) sin(psi_arg) 0; -sin(psi_arg) -cos(psi_arg) 0; 0 0 1];

DEBUG = false;

% Get P_i
[index_P_i] = get_P_i(Trajectoire, P_kite, last_index);
new_index = index_P_i;

% V_dir points towards the target
% V_dir = Trajectoire(index_P_i).position - P_kite;
% k0_V_dir = WR_M_k0(phi, theta).'*V_dir;


if(index_P_i == length(Trajectoire))
    index_next_P_i = 1;
else
    index_next_P_i = index_P_i+ 1;
end

% V_dir_2 points towards the target
V_dir_2 = Trajectoire(index_next_P_i).position - P_kite;
k0_V_dir_2 = WR_M_k0(phi, theta).'*V_dir_2;

% V_next points from the target P_i towards the next point P_i+1
V_next = Trajectoire(index_next_P_i).position - Trajectoire(index_P_i).position;
k0_V_next = WR_M_k0(phi, theta).'*V_next;


% V_psi is the current yaw (same as vector x_b)
V_psi = WR_M_k0(phi, theta)*k0_M_b(psi)*[8;0;0];
k0_V_psi = WR_M_k0(phi, theta).'*V_psi;

alpha = signed_angle_between_vectors(k0_V_next, k0_V_psi, [0;0;1]);
alpha_0 = signed_angle_between_vectors(k0_V_next, k0_V_dir_2, [0;0;1]);
alpha_V_dir = signed_angle_between_vectors([1;0;0], k0_V_next, [0;0;1]);

% Distance from the kite to the target P_i
error = norm(V_dir_2);
alpha_des = atan(error/d0)*alpha_0/(pi/2);

if(DEBUG)
    plot_vector(P_kite, V_dir, 'k', 'V_{dir}')
    plot_vector(P_kite, V_next, 'b', 'V_{next}')
    plot_vector(P_kite, V_psi, 'r', '\psi')
    disp(['alpha = ' num2str(floor(alpha*180/pi)) '°'])
    disp(['alpha_0 = ' num2str(floor(alpha_0*180/pi)) '°'])
    disp(['alpha_des = ' num2str(floor(alpha_des*180/pi)) '°'])
    plot_vector(P_kite, WR_M_k0(phi,theta)*[5;0;0],'m','x_{k0}')
    plot_vector(P_kite, WR_M_k0(phi,theta)*[0;5;0],'m','y_{k0}')
    plot_vector(P_kite, WR_M_k0(phi,theta)*[0;0;-5],'m','z_{k0}')
    V_des =  WR_M_k0(phi, theta)*[cos(alpha_V_dir+alpha_des); sin(alpha_V_dir+alpha_des); 0]*10;
    plot_vector(P_kite,V_des, '--r', 'V_{des}')
    plot3(Trajectoire(index_P_i).position(1), Trajectoire(index_P_i).position(2), Trajectoire(index_P_i).position(3), 'g.','markerSize',10)
end

end

