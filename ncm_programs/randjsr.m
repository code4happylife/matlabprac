function r = randjsr(p,q)
%RANDJSR   Shift register random number floating point generator.
%  Based on the parameters used by MATLAB version 4.
%  The statement
%     r = randjsr
%  generates a single uniformly distributed random number.
%  The statement
%     r = randjsr(m,n)
%  generates an m-by-n random matrix.
%  The statement
%     clear randjsr
%  will cause the generator to reinitialize itself.
%  The function can not accept any other starting seed.
%
%  See also RANDGUI, RANDSSP.

persistent j
if isempty(j)
   j = 123456789;
end

if nargin < 1, p = 1; end
if nargin < 2, q = p; end
r = zeros(p,q);
for k = 1:p*q
   j = bitxor(j,bitshift(j,13,32));
   j = bitxor(j,bitshift(j,-17,32));
   j = bitxor(j,bitshift(j,5,32));
   r(k) = pow2(j,-32);
end
