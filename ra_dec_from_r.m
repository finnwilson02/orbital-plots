function [ra dec] = ra_dec_from_r(r)

a = r(1)/norm(r);   %find the unit position vector for r
b = r(2)/norm(r);
c = r(3)/norm(r);

dec = asind(c);     %find the declination from the direction z

if b > 0            %find right ascension from the unit position vector
    ra = acosd(a/cosd(dec));    
else
    ra = 360 - acosd(a/cosd(dec));
end

end