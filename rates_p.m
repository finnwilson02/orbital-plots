function dydt = rates_p(t,f)

mu = 398.58e3;  %gravitational parameter (km^3 s^-2)

%take x,y,z components of position and velocity
x = f(1);
y = f(2);
z = f(3);
vx = f(4);
vy = f(5);
vz = f(6);

%calculate the x,y,z components of acceleration
r = norm([x y z]);
p = pertubation(f);
ax = -mu*x/r^3 + p(1);
ay = -mu*y/r^3 + p(2);
az = -mu*z/r^3 + p(3);

%return dydt
dydt = [vx vy vz ax ay az]';

end

