function Q = splinequad(x,y)
%SPLINEQUAD  Discrete quadrature with shape-preserving piecewise cubic
%  Q = splinequad(x,y) approximates int(f(x),x(1),x(n)) from the
%  data y(j) = f(x(j)), j = 1:n.
%
%  See SPLINE, PCHIPQUAD.
 
%  First derivatives
 
   h = diff(x);
   delta = diff(y)./h;
   d = splineslopes(h,delta);

%  Trapezoid rule, with spline correction.
   
   k = 1:length(h);
   T = sum(h.*(y(k+1)+y(k)))/2;
   D = sum(h.^2.*(d(k+1)-d(k)))/12;
   Q = T - D;

% -------------------------------------------------------

function d = splineslopes(h,delta)
%  SPLINESLOPES  Slopes for cubic spline interpolation.
%  splineslopes(h,delta) computes d(k) = S'(x(k)).
%  Uses not-a-knot end conditions.

%  Diagonals of tridiagonal system

   n = length(h)+1;
   a = zeros(size(h)); b = a; c = a; r = a;
   a(1:n-2) = h(2:n-1);
   a(n-1) = h(n-2)+h(n-1);
   b(1) = h(2);
   b(2:n-1) = 2*(h(2:n-1)+h(1:n-2));
   b(n) = h(n-2);
   c(1) = h(1)+h(2);
   c(2:n-1) = h(1:n-2);

%  Right-hand side

   r(1) = ((h(1)+2*c(1))*h(2)*delta(1)+h(1)^2*delta(2))/c(1);
   r(2:n-1) = 3*(h(2:n-1).*delta(1:n-2)+h(1:n-2).*delta(2:n-1));
   r(n) = (h(n-1)^2*delta(n-2)+(2*a(n-1)+h(n-1))*h(n-2)*delta(n-1))/a(n-1);

%  Solve tridiagonal linear system

   d = tridisolve(a,b,c,r);
