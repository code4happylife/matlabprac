function Q = pchipquad(x,y)
%PCHIPQUAD  Discrete quadrature with shape-preserving piecewise cubic
%  Q = pchipquad(x,y) approximates int(f(x),x(1),x(n)) from the
%  data y(j) = f(x(j)), j = 1:n.
%
%  See PCHIP, SPLINEQUAD.
 
%  First derivatives
 
   h = diff(x);
   delta = diff(y)./h;
   d = pchipslopes(h,delta);

%  Trapezoid rule, with pchip correction.
   
   k = 1:length(h);
   T = sum(h.*(y(k+1)+y(k)))/2;
   D = sum(h.^2.*(d(k+1)-d(k)))/12;
   Q = T - D;

% -------------------------------------------------------

function d = pchipslopes(h,delta)
%  PCHIPSLOPES  Slopes for shape-preserving Hermite cubic
%  interpolation.  pchipslopes(h,delta) computes d(k) = P'(x(k)).

%  Slopes at interior points
%  delta = diff(y)./diff(x).
%  d(k) = 0 if delta(k-1) and delta(k) have opposites signs
%         or either is zero.
%  d(k) = weighted harmonic mean of delta(k-1) and delta(k)
%         if they have the same sign.

   n = length(h)+1;
   d = zeros(size(h));
   k = find(sign(delta(1:n-2)).*sign(delta(2:n-1)) > 0) + 1;
   w1 = 2*h(k)+h(k-1);
   w2 = h(k)+2*h(k-1);
   d(k) = (w1+w2)./(w1./delta(k-1) + w2./delta(k));

%  Slopes at endpoints

   d(1) = pchipendpoint(h(1),h(2),delta(1),delta(2));
   d(n) = pchipendpoint(h(n-1),h(n-2),delta(n-1),delta(n-2));

% -------------------------------------------------------

function d = pchipendpoint(h1,h2,del1,del2)
%  Noncentered, shape-preserving, three-point formula.
   d = ((2*h1+h2)*del1 - h1*del2)/(h1+h2);
   if sign(d) ~= sign(del1)
      d = 0;
   elseif (sign(del1) ~= sign(del2)) & (abs(d) > abs(3*del1))
      d = 3*del1;
   end
