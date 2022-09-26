function [r,v] = coe2sv(coe,mu)

h = coe(1);    %angular momentum (kg^2/s)
e = coe(2);    %eccentricity
RA = coe(3);   %right ascension of the ascending node (rad)
incl = coe(4); %inclination of the orbit (rad)
w = coe(5);    %argument of perigee (rad)
TA = coe(6);   %true anomaly (rad)

%find position and velocity in the perifocal frame
rp = (h^2/mu) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] + sin(TA)*[0;1;0]);
vp = (mu/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);

%find the rotation matrices
%rotation around the z-axis through angle of right ascension
R3_W = [cos(RA) sin(RA) 0; -sin(RA) cos(RA) 0; 0 0 1];      
%rotation about x-axis though angle of inclination
R1_i = [1 0 0; 0 cos(incl) sin(incl); 0 -sin(incl) cos(incl)];
%rotation about z-axis through argument of perigee
R3_w = [cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1];

%use rotation matrices to find matrix of transformation from perifocal to
%geocentric equatorial frame
Q_pX = (R3_w*R1_i*R3_W)';

%use Q to find position and velocity in the geocentric equatorial frma
r = Q_pX*rp;
v = Q_pX*vp;

%make r,v row vectors
r = r';
v = v'

end