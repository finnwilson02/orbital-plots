function p = pertubation(f)

J2 = 0.00108263;
Re = 6371.1;    %radius of Earth (km)
mu = 398.58e3;  %gravitational parameter (km^3 s^-2)

x = f(1);
y = f(2);
z = f(3);
r = norm([x y z]);

a = (3/2)*(J2*mu*Re^2)/r^4;

i = (x/r)*(5*z^2/r^2 - 1);
j = (y/r)*(5*z^2/r^2 - 1);
k = (z/r)*(5*z^2/r^2 - 3);

p = a*[i,j,k];

end

