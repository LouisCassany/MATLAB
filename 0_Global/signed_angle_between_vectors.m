function signed_angle = signed_angle_between_vectors(v1, v2, n)
% Compute the angle (rad) between two vectors.
%
% signed_angle = signed_angle_between_vectors(v1, v2, n)
%
% v1, v2 are the 3D vectors between which the signed angle is to be 
% computed n is a vector used to sign the angle (right hand rule) 
% and can't be in the plane (v1,v2).
    
% Compute the dot product and cross product of the two vectors
dot_product = dot(v1, v2);
cross_product = cross(v1, v2);

% Compute the magnitude of the cross product
mag_cross_product = norm(cross_product);

% Compute the unsigned angle between the vectors using the dot product
unsigned_angle = atan2(mag_cross_product, dot_product);

% Compute the signed angle by checking the orientation of the reference
% vector
if dot(cross_product, n) < 0
    signed_angle = -unsigned_angle;
else
    signed_angle = unsigned_angle;
end
end

