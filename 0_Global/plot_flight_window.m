function plot_flight_window(r)
% Plot a flight window with a tether length r.

% quiver3(0,0,0,r,0,0,'r','HandleVisibility','off'); text(r,0,0,'X');
% hold on
% quiver3(0,0,0,0,r,0,'g','HandleVisibility','off'); text(0,r,0,'Y');
% quiver3(0,0,0,0,0,1); text(0,0,1,'z');

t=linspace(0,pi,20);
arc = [zeros(length(t),1) cos(t).' -sin(t).'];
My = @(alpha) [cos(alpha) 0 sin(alpha); 0 1 0; -sin(alpha) 0 cos(alpha)];
for theta = linspace(0,-pi/2,10)
    arc2 = My(theta)*arc.'*r;
    arc2 = arc2.';
    plot3(arc2(:,1), arc2(:,2), arc2(:,3),':k','HandleVisibility','off')
    hold on
end

t=linspace(0,pi/2,20);
arc = [zeros(length(t),1) cos(t).' -sin(t).'];
Mz = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
for theta = linspace(0,-pi,19)
    arc2 = Mz(theta)*arc.'*r;
    arc2 = arc2.';
    plot3(arc2(:,1), arc2(:,2), arc2(:,3),':k','HandleVisibility','off')
    hold on
end


hold off
a = gca;
a.XDir = 'Reverse';
a.ZDir = 'Reverse';
end

