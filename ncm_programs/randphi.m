function r = randphi(p,q)
% RANDPHI.  Another textbook random number generator.
% The statement
%    r = randphi
% generates a single random number.
% The statement
%    r = randphi(m,n)
% generates an m-by-n random matrix.
%
% See also RANDGUI, RANDMCG.

persistent x phi
if isempty(x)
   phi = (1+sqrt(5))/2;
   x = 0;
end

if nargin < 1, p = 1; end
if nargin < 2, q = p; end
r = zeros(p,q);
for k = 1:p*q
   x = rem(x + phi, 1);
   r(k) = x;
end
