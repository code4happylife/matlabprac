function [v,d] = perpchip(x,y,u)
%PERPCHIP  Periodic piecewise cubic Hermite interpolation.
%  v = perpchip(x,y,u) finds the shape preserving piecewise
%  cubic interpolant P(x), with P(x(j)) = y(j), and returns
%  v(k) = P(u(k)).
%
%  See PCHIPTX, SPLINETX.
 
%  First derivatives
 
   h = diff(x);
   delta = diff(y)./h;
   d = pchipslopes(h,delta);

%  Piecewise polynomial coefficients

   n = length(x);
   k = 1:n-1;
   c = (3*delta - 2*d(k) - d(k+1))./h;
   b = (d(k) - 2*delta + d(k+1))./h.^2;

%  Find subinterval indices, x(k) <= u < u(k+1)

   k = ones(size(u));
   for j = 2:n-1
      k(u >= x(j)) = j;
   end

%  Evaluate interpolant

   s = u - x(k);
   v = y(k) + s.*(d(k) + s.*(c(k) + s.*b(k)));


% -------------------------------------------------------

function d = pchipslopes(h,delta);
% PCHIPSLOPES  Slopes for periodic hermite cubic interpolation.
% pchipslopes(h,delta) computes d(k) = P'(x(k)).

%  Extend by periodicity

   n = length(h)+1;
   h(n) = h(1);
   delta(n) = delta(1);

%  Slopes at all points
%  delta = diff(y)./diff(x).
%  d(k) = 0 if delta(k-1) and delta(k) have opposites signs
%         or either is zero.
%  d(k) = weighted harmonic mean of delta(k-1) and delta(k)
%         if they have the same sign.

   d = zeros(size(delta));
   k = find(sign(delta(1:n-1)).*sign(delta(2:n)) > 0) + 1;
   w1 = 2*h(k)+h(k-1);
   w2 = h(k)+2*h(k-1);
   d(k) = (w1+w2)./(w1./delta(k-1) + w2./delta(k));
   d(1) = d(n);
