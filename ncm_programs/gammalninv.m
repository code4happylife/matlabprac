function y = gammalninv(x);
% GAMMALNINV   Inverse of log gamma function.
%   y = gammalninv(x) computes y so that gammaln(y) = x.
%   x must be > x0, then
%   y will be > y0, where
%   y0 = solve('Psi(x)') = 1.4616321449683623412626595423257
%   x0 = log(gamma(y0))  = -.12148629053584960809551455717769

x0 = double(-.12148629053584960809551455717769);
y0 = double(1.4616321449683623412626595423257);
y = zeros(size(x));
for k = 1:length(x(:))
   if x(k) < x0
      warning(['Function undefined for x < ' num2str(x0)])
      y(k) = NaN;
   else
      y(k) = fzerotx(@(y) gammaln(y)-x(k),[y0,x(k)/2+4]);
   end
end
