function [index_P_i, M_i] = get_P_i(Trajectoire, P_kite, index_last_P_i)
    % Réarrange la trajectoire afin de la parcourir en partant du point
    % d'index index_last_P_i
    Trajectoire_2 = [Trajectoire(index_last_P_i:end) Trajectoire(1:index_last_P_i)].';
    
    best_distance = inf;
  
    % Returned optimal point index
    index_P_i = index_last_P_i;
    
    % Check the next 5 points
    for index = 1 : 5
        next_P = Trajectoire_2(index).position;
        distance = norm(next_P - P_kite);
        if(distance < best_distance)
            best_distance = distance;
            index_P_i =  Trajectoire_2(index).index;
        end
    end
    
    % Get the transformation matrix
    [x_i, y_i, z_i] = plot_reperes(index_P_i, Trajectoire, 5, P_kite, false);
    
    M_i = [x_i y_i z_i];
end
