function finitesierpinski(n)
%FINITESIERPINSKI   Sierpinski's triangle.
%   FINITESIERPINSKI(N) does not run forever, but plots N points.
%   FINITEFERN with no arguments plots 5000 points.
%   See also FERN.

if nargin < 1, n = 100000; end

A = [1/2 0; 0 1/2];
b1 = [0 0]';
b2 = [1/2 0]';
b3 = [1/4 sqrt(3)/4]';

x = [0 0]';
xs = zeros(2,n);
xs(:,1) = x;
for j = 2:n
   r = ceil(3*rand);
   switch r
      case 1, x = A*x + b1;
      case 2, x = A*x + b2;
      case 3, x = A*x + b3;
   end
   xs(:,j) = x;
end

shg
set(gcf,'color','white')
darkgreen = [0 0 0];
if n <= 5000, dotsize = 6; else dotsize = 1; end
plot(xs(1,:),xs(2,:),'.','markersize',dotsize,'color',darkgreen);
axis([0 1 0 1])
axis square
axis off
