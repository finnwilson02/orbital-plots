function [r,v] = lagrange(r0,v0,theta)

mu = 398600;

%calculate magnitude of r0 and v0
r_vec = r0;
v_vec = v0;
r0 = sqrt(dot(r0,r0));
v0 = sqrt(dot(v0,v0));

%calculate the radial component of v0
vr = dot(v_vec,r_vec)/r0;

%calculate angular momentum h
h = r0*sqrt(v0^2 - vr^2);

%calculate r
temp = 1 + (h^2/(mu*r0)-1)*cos(theta) - ((h*vr)/mu)*sin(theta);
r = (h^2/mu)*1/temp;

%calulate f, f', g, g'
f = 1 - (mu*r)*(1-cos(theta))/h^2;
g = (r*r0)*sin(theta)/h;
fdot = (mu/h)*((1-cos(theta))/sin(theta))*(mu*(1-cos(theta))/h^2 - 1/r0 - 1/r);
gdot = 1 - mu*r0*(1-cos(theta))/h^2;

%calculate r and v
r = f*r_vec + g*v_vec;
v = fdot*r_vec + gdot*v_vec;

end