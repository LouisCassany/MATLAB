function [alpha, alpha_des, new_index] = asservissement_vitesse(P_kite, V_kite, Trajectoire, d0, last_index)
    [index_P_i, M_i] = get_P_i(Trajectoire, P_kite, last_index);
    new_index = index_P_i;
    V_dir = Trajectoire(index_P_i).position - P_kite;
    R_i_ref = M_i.'*M_i(:,3);
    R_i_x_i = M_i.'*M_i(:,1);
    R_i_V_dir = M_i.'*V_dir;
    R_i_V_dir(3) = 0;
    alpha_0 = signed_angle_between_vectors(R_i_x_i,R_i_V_dir,R_i_ref);
    R_i_V_kite = M_i.'*V_kite;
    R_i_V_kite(3) = 0;
    alpha = signed_angle_between_vectors(R_i_x_i,R_i_V_kite,R_i_ref);
    d = norm(V_dir);
    alpha_des = atan(d/d0)*alpha_0/(pi/2);
end