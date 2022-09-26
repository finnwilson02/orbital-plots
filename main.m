clear all;
clc;
%%
%Question 1
%Parameters and Constants
deg = pi/180;   %convert degrees to radians
hours = 3600;   %convert hours to seconds
mu = 398.58e3;  %gravitational parameter (km^3 s^-2)
Re = 6371.1;    %radius of Earth (km)

%Classical Orbital Elements
h = 52688.82;   %angular momentum (kg^2/s)
e = 0.0138867;  %eccentricity
RA = 328.6010;  %right ascension of the ascending node (rad)
incl = 78.3205; %inclination of the orbit (rad)
w = 40.3797;    %argument of perigee (rad)
TA = 320.2653;  %true anomaly (rad)

%Get State Vectors
coe = [h, e, RA*deg, incl*deg, w*deg, TA*deg];  %concatenate classical orbital elements
[r0 v0] = coe2sv(coe,mu);                       %call coe2sv to get state vectors expressing orbit

%Set up timespan and call ODE function
t0 = 0;
tf = 5786.62;
y0 = [r0 v0]';
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,y] = ode45(@rates,[t0 tf],y0,options);

%plot the orbit found with ode45
figure(10);
plotOrbit(t,tf,y,y0,'r');

%%
%Question 2
%Parameters and Constants
deg = pi/180;   %convert degrees to radians
hours = 3600;   %convert hours to seconds
mu = 398.58e3;  %gravitational parameter (km^3 s^-2)
Re = 6371.1;    %radius of Earth (km)
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

%Classical Orbital Elements
h = 52688.82;   %angular momentum (kg^2/s)
e = 0.0138867;  %eccentricity
RA = 328.6010;  %right ascension of the ascending node (rad)
incl = 78.3205; %inclination of the orbit (rad)
w = 40.3797;    %argument of perigee (rad)
TA = 0;  %true anomaly (rad)

%Get State Vectors
coe1 = [h, e, RA*deg, incl*deg, w*deg, TA*deg]; %concatenate classical orbital elements
[r0 v0] = coe2sv(coe1,mu);                      %call coe2sv to get state vectors expressing orbit

%Set up timespan and call ODE function
t0 = 0;             %starts in orbit A
tA = 5786.62;      %enters transfer orbit
tB = tA + 6520.01435*0.5;  %enters orbit B
tf = tB + 3286.29*2;  %end simulation at return to perigee
y0 = [r0 v0]';
[t1,y1] = ode45(@rates_p,[t0 tA],y0,options);

%Complete transfer orbit
rA = y1(end,1:3);
vA = y1(end,4:6);
v = norm([vA(1),vA(2),vA(3)]);
direction = vA./v;
dV_a = direction*0.26768;
vA = vA + dV_a;
yA = [rA vA]';
[t2,y2] = ode45(@rates_p,[tA tB],yA,options);

%Complete target orbit
h2 = 54813.5971;     %angular momentum (kg^2/s)
e2 = 0.077470;       %eccentricity
RA2 = 347.06;      %right ascension of the ascending node (rad)
incl2 = 100.4478;    %inclination of the orbit (rad)
w2 = 40.37197;        %argument of perigee (rad)
TA2 = 180;           %true anomaly (rad)

coe2 = [h2, e2, RA2*deg, incl2*deg, w2*deg, TA2*deg]; %concatenate classical orbital elements
[r2,v2] = coe2sv(coe2,mu);                      %call coe2sv to get state vectors expressing orbit
y0_2 = [r2,v2]';
[t3,y3] = ode45(@rates_p,[tB tf],y0_2,options);

%plot orbits on the same figure
figure(2);
plotOrbit(t1,tA,y1,y0,'#42b3a7');
hold on;
plotOrbit(t2,tB,y2,y0,'#ab3933');
hold on;
plotOrbit(t3,tf,y3,y0,'#6841b5');

%%
%Question 3
%initialise parameters
rP = 8171.1;
rA = 6996.1;
TA = 90;
W = 347.06;
i = 100.4478;
wp = 40.3797;
n = 5;

%call the ground track function
f = [rP,rA,TA,W,i,wp,n];
figure(3);
ground_track(f);

T = 1.82571*3600;       %period in seconds
days = 1.82571/24;      %convert to days
year = 365*24*3600;     %conversion from seconds to years
n1 = floor(year/T);     %create a time vector that goes up to the number
t = 1:n1;               %of orbital periods in a year
dRA_SSO = 0.0760658;    %pressecion of earth relative to sun
y = 328.601 + 0.075895*t;   %right ascension over time (earth fixed)
y1 = y - dRA_SSO*t;         %right ascension over time (sun fixed)

for i = 1:length(y)
    
    if y(i) >= 360          %keep right ascension between [0 360]
        y(i) = y(i) - 360;
    end
end

%plot right ascension relative to earth and the sun
figure(4);
plot(t*days,y);
yline(328.601,'-','Original Right Ascension');
ylabel('Right Ascension of Ascending Node (degrees)');
xlabel('Days Elapsed');

figure(5);
plot(t*days,y1);
ylabel('Right Ascension of Ascending Node (degrees)');
xlabel('Days Elapsed');
axis([0 360 0 400]);

%%
%Question 4
deg = pi/180;   %convert degrees to radians
hours = 3600;   %convert hours to seconds
mu = 398.58e3;  %gravitational parameter (km^3 s^-2)
Re = 6371.1;    %radius of Earth (km)

wgs84 = wgs84Ellipsoid('kilometer'); %create a reference ellipsoid
lat = -36.8509;               %initialise latitude, longitude and altitude
lon = 174.7645;
h = 42164-Re;

[X,Y,Z] = geodetic2ecef(wgs84,lat,lon,h);   %convert geodetic coordinates
r1 = [X,Y,Z];                               %to Earth-centered Earth-fixed
r2 = -r1 + 0.1;

T = 86164.06;
t = T/2;        %find the time between r1 and r2

[v1,v2] = lambert(r1,r2,t,'retro');   %use the lambert function to calculate
v1 = real(v1);                      %velocity at r1 and r2
coe = sv2coe(r1,v1,mu);             %use r1,v1 to find classic orbital elements

%display results to console
fprintf('–––––––––––––––––––––––––––––––––––––––––––––––––––––')
fprintf('\n Lambert''s Problem\n')
fprintf('\n\n Input data:\n');
fprintf('\n Gravitational parameter (km^3/s^2) = %g\n', mu);
fprintf('\n r1 (km) = [%g %g %g]',r1(1), r1(2), r1(3))
fprintf('\n r2 (km) = [%g %g %g]',r2(1), r2(2), r2(3))
fprintf('\n Elapsed time (s) = %g', t);
fprintf('\n\n Solution:\n')
fprintf('\n v1 (km/s) = [%g %g %g]',v1(1), v1(2), v1(3))
fprintf('\n v2 (km/s) = [%g %g %g]',v2(1), v2(2), v2(3))
fprintf('\n\n Orbital elements:')
fprintf('\n Angular momentum (km^2/s) = %g', coe(1))
fprintf('\n Eccentricity = %g', coe(2))
fprintf('\n Inclination (deg) = %g', coe(4)/deg)
fprintf('\n RA of ascending node (deg) = %g', coe(3)/deg)
fprintf('\n Argument of perigee (deg) = %g', coe(5)/deg)
fprintf('\n True anomaly initial (deg) = %g', coe(6)/deg)
fprintf('\n True anomaly final (deg) = %g', coe(6)/deg)
fprintf('\n Semimajor axis (km) = %g', coe(7))
fprintf('\n Periapse radius (km) = %g', coe(1)^2/mu/(1 + coe(2)))

%Set up timespan and call ODE function
t0 = 0;
t1 = 1.3184*86164.06;
t2 = t1 + 80526.3695;
t3 = t2 + 86164.06;
y0 = [r1 v1]';
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t_i,y1] = ode45(@rates_p,[t0 t1],y0,options);
figure(6);
plotOrbit(t_i,t1,y1,y0,'r');
hold on;

%enter classic orbital elements for phasing orbit
h_p = 126635.0036;
e_p = 0.041593;
RA_p = 0;
incl_p = 0;
w_p = (360-253.607)*deg;
TA_p = (180)*deg;
coe_p = [h_p,e_p,RA_p,incl_p,w_p,TA_p];

%convert to geocentric and plot phasing orbit
[r_p v_p] = coe2sv(coe_p,mu);
y_p = [r_p v_p]';
[t_p,y2] = ode45(@rates_p,[t1 t2],y_p,options);
plotOrbit(t_p,t2,y2,y_p,'b');

%enter orbital elements for the geostationary orbit
h_t = 129636.9049;      
e_t = 0;
RA_t = 0;
incl_t = 0;
w_t = (360-253.607)*deg;
TA_t = 180*deg;

%convert to geocentric, simulate one orbital period and plot results
coe_t = [h_t,e_t,RA_t,incl_t,w_t,TA_t];
[r_t v_t] = coe2sv(coe_t,mu);
y_t = [r_t v_t]';
[t_t,y3] = ode45(@rates_p,[t2 t3],y_t,options);
plotOrbit(t_t,t3,y3,y_t,'g');

%plot ground track for geosynchronous orbit
figure(7);
r_p1 = 42164;
r_a1 = 42164;
TA_1 = 174.136;
W_1 = 106.393;
incl_1 = 141.151;
wp_1 = 113.025;
n_1 = 1;
f1 = [r_p1,r_a1,TA_1,W_1,incl_1,wp_1,n_1];
ground_track(f1);

%plot ground track for geostationary orbit
figure(8);
r_p2 = 42164;
r_a2 = 42164;
TA_2 = 0;
W_2 = 0;
incl_2 = 0; 
wp_2 = 151.2093;
n_2 = 1;
f2 = [r_p2,r_a2,TA_2,W_2,incl_2,wp_2,n_2];
ground_track(f2);


