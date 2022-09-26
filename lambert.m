function [V1,V2] = lambert(R1,R2,t,str)

global r1 r2 A mu;
r1 = norm(R1);
r2 = norm(R2);
c12 = cross(R1, R2);
theta = acos(dot(R1,R2)/r1/r2);
mu = 398.58e3;  %gravitational parameter (km^3 s^-2)

if (~strcmp(str,'retro') && (~strcmp(str,'pro'))) %set default orbit to prograde
    string = 'pro';
    fprintf('\nPrograde trajectory assumed.\n')
end

if strcmp(str,'pro')             %keep TA between [0 2pi]
    if c12(3) <= 0
        theta = 2*pi - theta;
    end  
elseif strcmp(str,'retro')
    if c12(3) >= 0
    theta = 2*pi - theta;
    end
end


A = sin(theta)*sqrt(r1*r2/(1 - cos(theta)));    %calculate the constant, A

z = -100;
while F(z,t) < 0
    z = z + 0.1;        %find the starting value of z when F(z,t) changes signs
end

tol = 1.e-8;    %error tolerance
nmax = 5000;    %max iterations

ratio = 1;
n = 0;
while (abs(ratio) > tol) && (n <= nmax) %iterate until z is within error tolerance
    n = n + 1;
    ratio = F(z,t)/dFdz(z);
    z = z - ratio;
end

if n >= nmax        %error message if max iterations exceeded
    fprintf('\n\nNumber of iterations exceeds %g \n\n ',nmax)
end

f = 1 - y(z)/r1;        %calculate Lagrange coefficients
g = A*sqrt(y(z)/mu);
gdot = 1 - y(z)/r2;     %g' wrt time
V1 = 1/g*(R2 - f*R1);   %velocities at r1 and r2
V2 = 1/g*(gdot*R2 - R1);

end

function var = y(z)     %functions to simplify the lambert problem function
    global r1 r2 A mu;
    var = r1 + r2 + A*(z*S(z) - 1)/sqrt(C(z));
end

function var = F(z,t)
    global r1 r2 A mu;
    var = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mu)*t;
end

function var = dFdz(z)
    global r1 r2 A mu;
    if z == 0
        var = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
    else
        var = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) + A*sqrt(C(z)/y(z)));
    end
end

function var = C(z)
    var = stumpC(z);
end

function var = S(z)
    var = stumpS(z);
end
