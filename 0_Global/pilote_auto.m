function [control_signal, alpha, alpha_des, target_index, next_target_index, error] = pilote_auto(...
    theta_T0, phi_T0, n, a, Omega_T0, ...
    d0, P, r,...
    last_target_index, ...
    theta, phi, psi)

% Compute a current and desired orientation for the kite.
%
% [control_signal, alpha, alpha_des, target_index, next_target_index, error] = PILOTE_AUTO(theta_T0, phi_T0, n, a, Omega_T0, d0, P, r, last_target_index, theta, phi, psi)
%
% control_signal    : differential control (mm)
% alpha             : current orientation (rad)
% alpha_des         : desired orientation (rad)
% target_index      : target point index
% next_target_index : next target point index
% error             : Approx kite distance from the trajectory (m)
% ------------------------------------------------------------------
% theta_T0          : elevation of the trajectory center (rad)
% phi_T0            : azimuth of the trajectory center (rad)
% n                 : number of points for the trajectory
% a                 : trajectory width (m)
% Omega_T0          : trajectory rotation (rad)
% d0                : tracking aggresivity
% P                 : proportional controler
% r                 : tether length (m)
% last_target_index : last tracking point index
% theta             : kite elevation (rad)
% phi               : kite azimuth (rad)
% psi               : kite yaw (rad)


%% Build the trajectory
s = linspace(0, 2 * pi, n+1); % Curvilinear abscissa

% Trajectory center
Ch = [r * cos(theta_T0) * cos(phi_T0); r * cos(theta_T0) * sin(phi_T0); -r * sin(theta_T0)];

% Rotation matrices
My = @(alpha) [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
Mz = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

% Transformation from R_WR to R_k0 + rotation
WR_M_k0 = Mz(phi_T0) * My(theta_T0 - pi / 2) * Mz(Omega_T0);

yTdes = a * (sin(s)) ./ (1 + cos(s).^2);
xTdes = a * (sin(s) .* cos(s)) ./ (1 + cos(s).^2);

% Trajectory points in R_k_0
k0_Tdes = [xTdes; yTdes; zeros(1,length(xTdes))];

% Trajectory points in R_WR + offset of Ch
WR_Tdes = (WR_M_k0*k0_Tdes) + Ch;

% Pre allocation for speed and remove simulink warning
Trajectory = repmat(struct("position", [0;0;0], "index", 0), 1, length(k0_Tdes)-1);

for index = 1:length(WR_Tdes)-1
    Point = WR_Tdes(:, index);
    Trajectory(index).position = Point ./ norm(Point) * r;
    Trajectory(index).index = index;
end


% Trajectory is now a an array of structures of the form:
% [struct('position',[x;y;z],'index',1),...]

%% Compute the target points

% Kite position
P_kite = r * [(cos(phi) .* cos(theta)).'; (sin(phi) .* cos(theta)).'; -sin(theta).'];

% Reorder the trajectory points so that the array start with the point of
% index last_target_index
Trajecjory_new = [Trajectory(last_target_index:end) Trajectory(1:last_target_index)].';
Trajecjory_new(end) = [];

best_distance = inf;

target_index = last_target_index;

% Check the next 5 points
for index = 1 : 5
    next_P = Trajecjory_new(index).position;
    distance = norm(next_P - P_kite);
    if(distance < best_distance)
        best_distance = distance;
        % new_target_index is the index of the nearest forward point
        target_index =  Trajecjory_new(index).index;
    end
end

if(target_index == length(Trajectory))
    next_target_index = 1;
else
    next_target_index = target_index + 1;
end

%% Compute the desired yaw angle

% Transform matrices
WR_M_k0 = @(phi_arg, theta_arg) [cos(phi_arg)*sin(theta_arg) -sin(phi_arg) -cos(phi_arg)*cos(theta_arg); sin(phi_arg)*sin(theta_arg) cos(phi_arg) -sin(phi_arg)*cos(theta_arg); cos(theta_arg) 0 sin(theta_arg)];
k0_M_b = @(psi_arg) [-cos(psi_arg) sin(psi_arg) 0; -sin(psi_arg) -cos(psi_arg) 0; 0 0 1];

% V_dir points from the kite towards the next target
V_dir = Trajectory(next_target_index).position - P_kite;
k0_V_dir = WR_M_k0(phi, theta).'*V_dir;

% V_next points from the target towards the next target point
V_next = Trajectory(next_target_index).position - Trajectory(target_index).position;
k0_V_next = WR_M_k0(phi, theta).'*V_next;

% V_psi is the current yaw (same as vector x_b)
V_psi = WR_M_k0(phi, theta)*k0_M_b(psi)*[10;0;0];
k0_V_psi = WR_M_k0(phi, theta).'*V_psi;

% Current orientation
alpha = signed_angle_between_vectors(k0_V_next, k0_V_psi, [0;0;1]);
alpha_0 = signed_angle_between_vectors(k0_V_next, k0_V_dir, [0;0;1]);

% Approx distance from the trajectory
error = norm(V_dir);

% Desired orientation
alpha_des = atan(error/d0)*alpha_0/(pi/2);

control_signal = P*(alpha_des-alpha);