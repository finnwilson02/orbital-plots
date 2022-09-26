function coe = sv2coe(R,V,mu)

eps = 1.e-10;                       %min eccentricity without rounding

r = norm(R);                        %find scalars for position and velcity
v = norm(V);

vr = dot(R,V)/r;                    %radial velocity

H = cross(R,V);                     %angular momentum vector and scalar
h = norm(H);

incl = acos(H(3)/h);                %inclination

N = cross([0 0 1],H);               %node line
n = norm(N);

if n ~= 0                           %calculate right ascension of ascending node
    RA = acos(N(1)/n); 
    
    if N(2) < 0
        RA = 2*pi - RA;             %keep right ascension within [0 2pi]
    end
else
    RA = 0;
end

E = 1/mu*((v^2 - mu/r)*R - r*vr*V); %eccentricity
e = norm(E);

if n ~= 0                           %argument of perigree
    if e > eps
        w = acos(dot(N,E)/n/e);
        
        if E(3) < 0                 %keep between [0 2pi]
            w = 2*pi - w;
        end
    else
    w = 0;                          %for a circular orbit, argument of perigee is 0
    end
else
    w = 0;
end

if e > eps                          %true anomaly
    TA = acos(dot(E,R)/e/r);
    
    if vr < 0
        TA = 2*pi - TA;             %keep between [0 2pi]
    end
else
    cp = cross(N,R);                %true anomaly for circular orbit
    
    if cp(3) >= 0
        TA = acos(dot(N,R)/n/r);
    else
        TA = 2*pi - acos(dot(N,R)/n/r);
    end
end

a = h^2/mu/(1 - e^2);               %semimajor axis

coe = [h e RA incl w TA a];         %concatenate results

end