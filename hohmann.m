%On: h, e, RA, incl, w, TA
function hohmann(O1,O2)

mu = 398600;

%find the perigee of orbit 1 and apogee of orbit 2
h1 = O1(1);
e1 = O1(2);
a1 = h1^2/(mu*(1-e1^2));
rp = a1 - a1*e1;
v1p = h1/rp;

h2 = O2(1);
e2 = O2(2);
a2 = h2^2/(mu*(1-e2^2));
ra = a2 + a2*e2;
v2a = h2/ra;

%find the angular momentum and velocity at apse points for orbit 3
h3 = sqrt((2*mu*rp*ra)/(ra + rp));
v3a = h3/ra;
v3p = h3/rp;

%calculate the delta v requirements
dv_p = v3p - v1p
dv_a = v2a - v3a

end

end