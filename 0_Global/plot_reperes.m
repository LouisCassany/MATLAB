function [x_i, y_i, z_i] = plot_reperes(index, Trajectoire, scale, Point, draw_plot)
    P_i = Trajectoire(index).position;

    % PREVENT BUG WHERE index_suivant HAPPENS TO BE 0
    if (index == length(Trajectoire))
        index_suivant = 1;
    else
        index_suivant = index + 1;
    end

    P_i_suivant = Trajectoire(index_suivant).position;
    x_i = (P_i_suivant - P_i) / norm(P_i_suivant - P_i);
    z_i = P_i_suivant / norm(P_i_suivant).';
    y_i = cross(z_i, x_i) / norm(cross(z_i, x_i)).';
    x_i = cross(y_i, z_i) / norm(cross(y_i, z_i));
    if(draw_plot)
        x_i = x_i * scale;
        y_i = y_i * scale;
        z_i = z_i * scale;
        plot_vector(Point, x_i, 'b');
        plot_vector(Point, y_i, 'r');
        plot_vector(Point, z_i, 'g');
    end
end