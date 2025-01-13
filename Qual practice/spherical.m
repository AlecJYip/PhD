function [r,theta,phi] = spherical(x,y,z)

r = sqrt(x^2+y^2+z^2);

rho = sqrt(x^2+y^2);

theta = atan(rho/z);

phi = atan2(y,x);


end

