% HANDINTERP

% Load data.

load myhand.dat
x = myhand(:,1);
y = myhand(:,2);
n = length(x);

% Interpolate using index as independent variable.

s = (1:n)';
t = (1:1/8:n)';
u = splinetx(s,x,t)
v = splinetx(s,y,t);

% Plot

figure(1)
clf
p = plot(x,y,'bo',u,v,'g-');
set(p(1),'markersize',4)
set(p(2),'color',[0 2/3 0])
axis([-.05 1.05 -.05 1.05])

% Interpolate in polar coordinates.
% Translate to make starlike with respect to origin.

figure(2)
x0 = (x(1) + x(n))/2;
y0 = (y(1) + y(n))/2;
x = x - x0;
y = y - y0;
w = w - (x0 + i*y0);
theta = atan2(y,x);
r = sqrt(x.^2 + y.^2);
plot(-theta,r,'.-')

figure(3)
polar(theta,r,'.-')

figure(4)
t = (theta(1):(theta(n)-theta(1))/256:theta(n))';
Z = pchiptx(theta,r,t).*exp(i*t);
S = splinetx(theta,r,t).*exp(i*t);
p = plot(x,y,'bo', ...
   real(w),imag(w),'g-', ...
   real(Z),imag(Z),'m-', ...
   real(S),imag(S),'k-');
set(p(1),'markersize',4)
set(p(2),'color',[0 2/3 0])
axis(ax)
