function plot_vector(P, V, style, textString)
% Wrapper for quiver3.
%
% plot_vector(P,V,style, optional : textString)
%
% P             : point to plot the vector
% V             : vector
% style         : style e.g 'k--'
% textString    : optionan text to display


if(length(P) == 3)
    %3D case
    if(nargin == 4)
        text(P(1)+V(1),P(2)+V(2),P(3)+V(3),textString,'Interpreter','tex');
    end
    quiver3(P(1), P(2), P(3), V(1), V(2), V(3), style,'AutoScale','off');
else
     %3D case
    if(nargin == 4)
        text(P(1)+V(1),P(2)+V(2),0,textString,'Interpreter','tex');
    end
    quiver3(P(1), P(2), 0, V(1), V(2), 0, style,'AutoScale','off');
end
  
end

