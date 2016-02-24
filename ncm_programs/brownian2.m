function brownian2(N,p)
% BROWNIAN2  2-D Brownian motion.
%   BROWNIAN(N,p) takes n steps with p points.
clf
shg
set(gcf,'doublebuffer','on')
delta = .002;
if nargin < 1, N = 10000; end
if nargin < 2, p = 100; end
x = zeros(p,1);
y = zeros(p,1);
h = plot(x,y,'.','erasemode','xor');
axis([-1 1 -1 1])
axis square
r = zeros(N,1);
R = zeros(N,1);
for n = 1:N
   x = x + delta*randn(p,1);
   y = y + delta*randn(p,1);
   set(h,'xdata',x,'ydata',y)
   drawnow
   r(n,1) = mean(sqrt(x.^2 + y.^2));
   R(n,1) = max(sqrt(x.^2 + y.^2));
end
pause
clf
n = (1:N)';
c = sqrt(n)\r;
C = sqrt(n)\R;
plot(n,[c*sqrt(n) C*sqrt(n) r R])
