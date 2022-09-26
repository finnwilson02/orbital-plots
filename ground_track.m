function ground_track(f);

%...Constants
deg = pi/180;
mu = 398580;  %gravitational parameter (km^3 s^-2)
J2 = 0.00108263;
Re = 6378;
we = (2*pi + 2*pi/365.26)/(24*3600);

%parse data from input
rP = f(1);          %perigee altitude (km)
rA = f(2);          %apogee altitude (km)
TAo = f(3)*deg;     %true anomaly (radians)
Wo = f(4)*deg;      %right ascension of ascending node (radians)
incl = f(5)*deg;    %inclination (radians)
wpo = f(6)*deg;     %argument of perigee (radians)
n_periods = f(7);   %number of periods

%compute classical orbital elements, time, and rates of node
%regression/perigee advance
a = (rA + rP)/2;
T = 2*pi/sqrt(mu)*a^(3/2);
e = (rA - rP)/(rA + rP);
h = sqrt(mu*a*(1 - e^2));
E = 2*atan(tan(TAo/2)*sqrt((1-e)/(1+e)));
M = E - e*sin(E);
to = M*(T/2/pi);
tf = to + n_periods*T;
fac = -3/2*sqrt(mu)*J2*Re^2/(1-e^2)^2/a^(7/2);
Wdot = fac*cos(incl);
wpdot = fac*(5/2*sin(incl)^2 - 2);

%propogate the orbit over the given time, converting to an earth-fixed
%frame of reference and computing RA, dec over the time interval
times = linspace(to,tf,1000);
ra = [];
dec = [];
theta = 0;

for i = 1:length(times)
    
    t = times(i);
    M = 2*pi/T*t;
    E = kepler_E(e, M);
    TA = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
    r = h^2/mu/(1 + e*cos(TA))*[cos(TA) sin(TA) 0]';
    W = Wo + Wdot*t;
    wp = wpo + wpdot*t;
    R1 = [cos(W) sin(W) 0; -sin(W) cos(W) 0; 0 0 1];
    R2 = [1 0 0; 0 cos(incl) sin(incl); 0 -sin(incl) cos(incl)];
    R3 = [cos(wp) sin(wp) 0; -sin(wp) cos(wp) 0; 0 0 1];
    QxX = (R3*R2*R1)';
    R = QxX*r;
    theta = we*(t - to);
    Q = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
    r_rel = Q*R;
    [alpha delta] = ra_dec_from_r(r_rel);
    ra = [ra; alpha];
    dec = [dec; delta];
end

%break up the ground track into seperate curves, each in the range [0 360]
%degrees
tol = 100;
curve_no = 1;
n_curves = 1;
k =0;
ra_prev = ra(1);

for i = 1:length(ra)
    if abs(ra(i) - ra_prev) > tol
        curve_no = curve_no + 1;
        n_curves = n_curves + 1;
        k = 0;
    end
    
    k = k + 1;
    RA{curve_no}(k) = ra(i);
    Dec{curve_no}(k) = dec(i);
    ra_prev = ra(i);
end

%plot the ground track
hold on
img = imread('map.jpg');
image('CData',img,'XData',[0 360],'YData',[90 -90]);
xlabel('East longitude (degrees)')
ylabel('Latitude (degrees)')
axis equal
grid on
for i = 1:n_curves
    plot(RA{i}, Dec{i},'LineWidth',2);
end
axis ([0 360 -90 90])
text( ra(1), dec(1), 'o Start')
text(ra(end), dec(end), 'o Finish')
line([min(ra) max(ra)],[0 0], 'Color','k')

%print orbital data
coe = [h e Wo incl wpo TAo];
[ro, vo] = coe2sv(coe, mu);
fprintf('\n ––––––––––––––––––––––––––––––––––––––––––––––––––––\n')
fprintf('\n Angular momentum = %g km^2/s' , h)
fprintf('\n Eccentricity = %g' , e)
fprintf('\n Semimajor axis = %g km' , a)
fprintf('\n Perigee radius = %g km' , rP)
fprintf('\n Apogee radius = %g km' , rA)
fprintf('\n Period = %g hours' , T/3600)
fprintf('\n Inclination = %g deg' , incl/deg)
fprintf('\n Initial true anomaly = %g deg' , TAo/deg)
fprintf('\n Time since perigee = %g hours' , to/3600)
fprintf('\n Initial RA = %g deg' , Wo/deg)
fprintf('\n RA_dot = %g deg/period' , Wdot/deg*T)
fprintf('\n Initial wp = %g deg' , wpo/deg)
fprintf('\n wp_dot = %g deg/period' , wpdot/deg*T)
fprintf('\n')
fprintf('\n r0 = [%12g, %12g, %12g] (km)', ro(1), ro(2), ro(3))
fprintf('\n magnitude = %g km\n', norm(ro))
fprintf('\n v0 = [%12g, %12g, %12g] (km)', vo(1), vo(2), vo(3))
fprintf('\n magnitude = %g km\n', norm(vo))
fprintf('\n ––––––––––––––––––––––––––––––––––––––––––––––––––––\n')

end