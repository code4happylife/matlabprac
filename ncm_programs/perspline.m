function v = perspline(x,y,u)
%PERSPLINE  Periodic spline function.
%  v = perspline(x,y,u) finds the piecewise cubic interpolatory
%  spline S(x), with S(x(j)) = y(j), and returns v(k) = S(u(k)).
%
%  See PERPCHIP, SPLINE, PCHIPTX.

%  First derivatives

   h = diff(x);
   delta = diff(y)./h;
   d = splineslopes(h,delta);

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

%  Evaluate spline

   s = u - x(k);
   v = y(k) + s.*(d(k) + s.*(c(k) + s.*b(k)));


% -------------------------------------------------------

function d = splineslopes(h,delta);
% SPLINESLOPES  Slopes for periodic cubic spline interpolation.
% splineslopes(h,delta) computes d(k) = S'(x(k)).

%  Set up sparse almost tridiagonal system and right hand side.

   rows = (size(h,2) > size(h,1));
   n = length(h)+1;
   h = [h(:); h(1)];
   delta = [delta(:); delta(1)];
   k = 1:n-1;
   A = spdiags([h(k+1) 2*(h(k+1)+h(k)) h(k)], -1:1, n-1, n-1);
   A(1,n-1) = h(n);
   A(n-1,1) = h(1);
   r = 3*(h(k+1).*delta(k)+h(k).*delta(k+1));

%  Solve almost tridiagonal linear system

   d = A\r;
   d = [d(n-1); d];
   if rows, d = d'; end
