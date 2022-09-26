function E = kepler_E(e,M)

error = 1.e-8; %set error tolerance

if M < pi           %set an initial value for E
    E = M + e/2;
else
    E = M - e/2;
end

ratio = 1;

while abs(ratio) > error    %iterate an approximation for E until it satisfies set tolerance
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;
end

end