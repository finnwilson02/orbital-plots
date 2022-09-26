function [ra dec] = find_ra_and_dec(tf)

times = linspace(to,tf,1000);
ra = [];
dec = [];
theta = 0;
for i = 1:length(times)
t = times(i);
M = 2*pi/T*t;
E = kepler_E(e, M);
TA = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
r =h^2/mu/(1 + e*cos(TA))*[cos(TA) sin(TA) 0]’;