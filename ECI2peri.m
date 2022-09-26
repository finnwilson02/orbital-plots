function PN = ECI2peri(w,RA,incl)

%find the rotation matrices
%rotation around the z-axis through angle of right ascension
R3_W = [cos(RA) sin(RA) 0; -sin(RA) cos(RA) 0; 0 0 1];      
%rotation about x-axis though angle of inclination
R1_i = [1 0 0; 0 cos(incl) sin(incl); 0 -sin(incl) cos(incl)];
%rotation about z-axis through argument of perigee
R3_w = [cos(w) sin(w) 0; -sin(w) cos(w) 0; 0 0 1];

PN = (R3_w*R1_i*R3_W);

end