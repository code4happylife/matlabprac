function brownian3(N,p)
% BROWNIAN3  3-D Brownian motion.
%   BROWNIAN3(n,p) takes n steps with p points.
clf
shg
set(gcf,'doublebuffer','on')
if nargin < 1, N = 10000; end
if nargin < 2, p = 100; end
delta = .005;
x = zeros(p,1);
y = zeros(p,1);
z = zeros(p,1);
h = plot3(x,y,z,'.','erasemode','xor','markersize',16);
box on
axis([-1 1 -1 1 -1 1])
daspect([1 1 1])
r = zeros(N,1);
R = zeros(N,1);
for n = 1:N
   x = x + delta*randn(p,1);
   y = y + delta*randn(p,1);
   z = z + delta*randn(p,1);
   set(h,'xdata',x,'ydata',y,'zdata',z)
   drawnow
   r(n,1) = mean(sqrt(x.^2 + y.^2 + z.^2));
   R(n,1) = max(sqrt(x.^2 + y.^2 + z.^2));
end
pause
clf
n = (1:N)';
c = sqrt(n)\r;
C = sqrt(n)\R;
plot(n,[c*sqrt(n) C*sqrt(n) r R])
