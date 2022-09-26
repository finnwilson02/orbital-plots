function plotOrbit(t,tf,y,y0,colour)

%find r (distance from center of Earth to satellite) at each t
for i = 1:length(t)
    r(i) = norm([y(i,1) y(i,2) y(i,3)]);
end

%find min and max distance and velocity at those points
[rmax imax] = max(r);
[rmin imin] = min(r);
v_rmax = norm([y(imax,4) y(imax,5) y(imax,6)]);
v_rmin = norm([y(imin,4) y(imin,5) y(imin,6)]);

r0 = y0(1:3);
v0 = y0(4:6);
hours = 3600;   %convert hours to seconds
Re = 6378;    %radius of Earth (km)

%Output to the command window:
fprintf('\n\n--------------------------------------------------------\n');
fprintf('\n Earth Orbit\n');
fprintf(' %s\n', datestr(now));
fprintf('\n The initial position is [%g, %g, %g] (km).',r0(1), r0(2), r0(3));
fprintf('\n Magnitude = %g km\n', norm(r0));
fprintf('\n The initial velocity is [%g, %g, %g] (km/s).',v0(1), v0(2), v0(3));
fprintf('\n Magnitude = %g km/s\n', norm(v0));
fprintf('\n Initial time = %g h.\n Final time = %g h.\n',0,tf/hours);
fprintf('\n The minimum altitude is %g km at time = %g h.',rmin-Re, t(imin)/hours);
fprintf('\n The speed at that point is %g km/s.\n', v_rmin);
fprintf('The velocity vector at that point is [%g %g %g]\n',y(imin,4),y(imin,5),y(imin,6));
fprintf('The location vector at that point is [%g %g %g]\n',y(imin,1),y(imin,2),y(imin,3));
fprintf('\n The maximum altitude is %g km at time = %g h.',rmax-Re, t(imax)/hours);
fprintf('\n The speed at that point is %g km/s\n', v_rmax);
fprintf('The velocity vector at that point is [%g %g %g]\n',y(imax,4),y(imax,5),y(imax,6));
fprintf('The velocity vector at that point is [%g %g %g]\n',y(imax,1),y(imax,2),y(imax,3));
fprintf('\n---------------------------------------------------------\n\n');


%Draw the Planet
earth_sphere('km'); %credit Will Campbell (2021)

%Draw the Axis
line([0 2*Re], [0 0], [0 0]); text(2*Re, 0, 0, 'X');
line( [0 0], [0 2*Re], [0 0]); text( 0, 2*Re, 0, 'Y');
line( [0 0], [0 0], [0 2*Re]); text( 0, 0, 2*Re, 'Z');

%Plot the Orbit and Label Initial/Final Points
hold on
plot3( y(:,1), y(:,2), y(:,3),'color',colour,'LineWidth',2);
line([0 r0(1)], [0 r0(2)], [0 r0(3)]);
plot3( y(1,1), y(1,2), y(1,3), 'o','color',colour,'LineWidth',2);
view([1,1,.4]);
xlabel('km');
ylabel('km');
zlabel('km');
end