function pendulum
% PENDULUM  Nonlinear, single pendulum.

L = 30;     % cm
g = 981;    % cm/s^2

Tsmall = 2*pi*sqrt(L/g)
disp(' ')
format short e
format compact

theta = [.001 (1:9)/10 (91:99)/100 (991:999)/1000]*pi;
T = zeros(size(theta));
Te = zeros(size(theta));
for k = 1:length(theta)
   T(k) = 4*sqrt(L/g)*myellipke(sin(theta(k)/2)^2);
   Te(k) = 4*sqrt(L/g)*ellipke(sin(theta(k)/2)^2);
   disp([theta(k),T(k),T(k)-Te(k)])
end
figure(1)
plot(theta,T,'-');
line([0 pi],[Tsmall Tsmall],'linestyle',':','color',[0 2/3 0])
axis([0 pi 0 5])
xlabel('Initial angle, \theta')
ylabel('Period, T')
title('Nonlinear pendulum')
pause

disp(' ')
theta = [.001 (1:7)/8 999/1000]'*pi;
tol = 1.e-6;
opts = odeset('reltol',tol,'abstol',tol);
figure(2)
c = [0*theta 2/3*(pi-theta)/pi 2/3*theta/pi];
for k = 1:length(theta)
   T = 4*sqrt(L/g)*ellipke(sin(theta(k)/2)^2);
   y0 = [theta(k) 0];
   [t,y] = ode45(@pendulode,[0 T],y0,opts,L,g);
   disp([theta(k),T,norm(y0-y(end,:))])
   plot(y(:,1),y(:,2),'-','color',c(k,:));
   hold on
   drawnow
end
hold off
title('Nonlinear pendulum')
xlabel('theta')
ylabel('thetadot')


% ------------------------

function ydot = pendulode(t,y,L,g);
ydot = [y(2); -(g/L)*sin(y(1))];


% ------------------------

function K = myellipke(ssquared)
K = quadtx(@ellipticintegrand,eps,1-eps,1.e-6,ssquared);


% ------------------------

function f = ellipticintegrand(t,ssquared)
f = 1./(sqrt(1-ssquared*t.^2).*sqrt(1-t.^2));

