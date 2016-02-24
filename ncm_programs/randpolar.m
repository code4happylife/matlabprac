function z = randpolar(p,q)
%RANDPOLAR   Polar algorithm for normally distributed random numbers.
%  The statement
%     r = randpolar
%  generates a single normally distributed random number.
%  The statement
%     r = randpolar(m,n)
%  generates an m-by-n normally distributed random matrix.
%  This function uses RAND, so resetting the state of RAND
%  will reset the state of RANDPOLAR.
%
%  See also RANDN, RAND.

if nargin < 1, p = 1; end
if nargin < 2, q = p; end
% Generate two at a time, so make length even.
pq = 2*ceil(p*q/2);
z = zeros(pq,1);
for k = 1:2:pq
   r = Inf;
   while r > 1
      u = 2*rand(2,1)-1;
      r = u'*u;
   end
   v = sqrt(-2*log(r)/r)*u;
   z(k:k+1) = v;
end
if pq > p*q, z(pq) = []; end
z = reshape(z,p,q);
