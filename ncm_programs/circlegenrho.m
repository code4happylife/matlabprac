function [rho, kappa] = circlegenrho(h)
if nargin < 1
   h = sqrt(2*(1-cos(2*pi/30)))
end
x = 1;
y = 0;
rmax = 1;
rmin = 1;
for n = 1:10000
   x = x + h*y;
   y = y - h*x;
   r = norm([x; y]);
   rmax = max(rmax,r);
   rmin = min(rmin,r);
end
rho = rmax/rmin;

A = [1 h; -h 1-h^2]
[X,D] = eig(A);
kappa = cond(X);
