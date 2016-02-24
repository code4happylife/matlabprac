function [tout,yout] = myrk4(F,tspan,y0,h,varargin)
%MYRK4  My own version of classical fourth order Runge-Kutta.
%   Usage is same as ODE23TX, except fourth argument is fixed step size, h.
%
%   MYRK4(F,TSPAN,Y0,H) with TSPAN = [T0 TFINAL] integrates the system
%   of differential equations y' = f(t,y) from t = T0 to t = TFINAL.
%   The initial condition is y(T0) = Y0.
%
%   With two output arguments, [T,Y] = MYRK4(...) returns a column 
%   vector T and an array Y where Y(:,k) is the solution at T(k).
%
%   With no output arguments, MYRK4 plots the emerging solution.
%
%   More than four input arguments, MYRK4(F,TSPAN,Y0,H,P1,P2,...),
%   are passed on to F, F(T,Y,P1,P2,...).

% Initalize variables.

t0 = tspan(1);
tfinal = tspan(2);
plotit = (nargout == 0);
t = t0;
y = y0(:);

% Initialize output.

if plotit
   odeplot(tspan,y,'init');
else
   tout = t;
   yout = y.';
end

% The main loop.

while t < tfinal

   if 1.1*abs(h) >= abs(tfinal - t)
      h = tfinal - t;
   end
  
   s1 = F( t, y, varargin{:});
   s2 = F( t+h/2, y+h/2*s1, varargin{:});
   s3 = F( t+h/2, y+h/2*s2, varargin{:});
   s4 = F( t+h, y+h*s3, varargin{:});
   t = t + h;
   y = y + h*(s1 + 2*s2 + 2*s3 + s4)/6;

   if plotit
      odeplot(t,y,'');
   else
      tout(end+1,1) = t;
      yout(end+1,:) = y.';
   end
end

if plotit
   odeplot([],[],'done');
end
