function Trajectoire = get_trajectory(theta_T0, phi_T0, n, a, Omega_0)
% Compute an eight shaped trajectory. The output is an array with the
% trajectory points.
%
% Trajectoire = GET_TRAJECTORY(theta_T0, phi_T0, n, a, Omega_0)
%
% Trajectoire   : [struct('position',[x;y;z],'index',1), ...] 
% theta_T0      : elevation of the trajectory center (rad)
% phi_T0        : azimuth of the trajectory center (rad)
% n             : number of points for the trajectory
% a             : arbitrary constant to control trajectory size
% Omega_0       : trajectory rotation (rad)

s = linspace(0, 2 * pi, n+1);

[r, ~] = init_model();

% Centre de la trajectoire, coordonnées polaires
Ch = [r * cos(theta_T0) * cos(phi_T0); r * cos(theta_T0) * sin(phi_T0); -r * sin(theta_T0)];

% Matrice de passages
Mx = @(alpha) [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];
My = @(alpha) [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
Mz = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

% Matrice de passage de R_WR -> R_k0
WR_M_k0 = Mz(phi_T0) * My(theta_T0 - pi / 2);

% Matrice de passage de R_WR -> R_k0 + rotation
WR_M_k0 = Mz(phi_T0) * My(theta_T0 - pi / 2) * Mz(Omega_0);

yTdes = a * (sin(s)) ./ (1 + cos(s).^2);
xTdes = a * (sin(s) .* cos(s)) ./ (1 + cos(s).^2);
k0_Tdes = [xTdes; yTdes; zeros(1, length(xTdes))];

% Pre allocation for speed and remove simulink warning
WR_Tdes = zeros(3,length(k0_Tdes));

for index = 1:length(k0_Tdes)
    P = k0_Tdes(:,index);
    WR_Tdes(:,index) = WR_M_k0 * P + Ch;
end

% Pre allocation for speed and remove simulink warning
Trajectoire = repmat(struct("position", [0;0;0], "index", 0), 1, length(k0_Tdes)-1);

for index = 1:length(WR_Tdes)-1
    Point = WR_Tdes(:, index);
    Trajectoire(index).position = Point ./ norm(Point) * r;
    Trajectoire(index).index = index;
end

end

