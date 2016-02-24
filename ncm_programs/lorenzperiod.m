function lorenzperioddriver
% Period of Lorenz attractor for several different values of rho.

rho = 99.65
p = lorenzperiod(rho,120,130)
pause

rho = 100.5
p = lorenzperiod(rho,80,90)
pause

rho = 160
p = lorenzperiod(rho,50,60)
pause

rho = 350
p = lorenzperiod(rho,40,44)


% ------------------------------------

function p = lorenzperiod(rho,t1,t2)
% lorenzperiod(rho,t1,t2).
% Integrate over tspan = [0 t2], but only display t = [t1 t2].
% t1 should be large enough for initial transients to die out.
% t2-t1 should be large enough to sample several periods.

sigma = 10;
beta = 8/3;
eta = sqrt(beta*(rho-1));
A = [ -beta    0     eta
         0  -sigma   sigma 
      -eta   rho    -1  ];

yc = [rho-1 eta eta];
y0 = yc + [0 0 3];
opts = odeset('events',@lorenzgstop,'reltol',1.e-8,'abstol',1.e-8);
[t,y,te,ye] = ode45(@lorenzeqn, [0 t2], y0, opts, A);
k = find(t>t1);
t = t(k);
y = y(k,:);
ke = find(te>t1);
te = te(ke);
ye = ye(ke,:);

% The event handling has found local maxima of norm(y(t)).
% A full period contains several of these local maxima.
% Identify maxima where all three components of y are nearly equal.
% The time between these distinguished maxima is the period.

ynorm = sqrt(sum(y.*y,2));
yenorm = sqrt(sum(ye.*ye,2));
m = max(find(yenorm > .999*max(yenorm)));
yem = ye(m,:);
n = length(te);
e = sum(abs(ye - yem(ones(n,1),:)),2);
k = e < 1.e-4*yenorm;
% disp([te ye e k])
d = diff(te(k));
if length(d) < 2 | (max(d)-min(d))/max(d) > .999
   warning('Has not reached steady state')
else
   minmeanmax = [min(d) mean(d) max(d)];
end

% The plots should not show any transient behavior.
% Dark green dots are local maxima of norm(y).
% Dark red dots show a complete period.

figure(1)
plot3(y(:,1),y(:,2),y(:,3),'-',ye(:,1),ye(:,2),ye(:,3),'.');
line(yem(1),yem(2),yem(3),'marker','.','markersize',18,'color',[2/3 0 0]);
line(yc(1),yc(2),yc(3),'marker','o','color',[0 2/3 0])
line(yc(1),-yc(2),-yc(3),'marker','o','color',[0 2/3 0])

figure(2)
for j = 1:3
   subplot(4,1,j)
   plot(t,y(:,j),'-',te,ye(:,j),'.',te(k),ye(k,j),'.')
end
subplot(4,1,4)
plot(t,ynorm,'-',te,yenorm,'.',te(k),yenorm(k),'.')
   
p = mean(d);

% ------------------------------------

function ydot = lorenzeqn(t,y,A)
%LORENZEQN  Equation of the Lorenz chaotic attractor.
%   ydot = lorenzeqn(t,y,A).
%   The differential equation is written in almost linear form.
%      ydot = A*y
%   where
%      A = [ -beta    0     y(2)
%               0  -sigma   sigma 
%            -y(2)   rho    -1  ];

A(1,3) = y(2);
A(3,1) = -y(2);
ydot = A*y;


% ------------------------------------

function [gstop,isterm,direct] = lorenzgstop(t,y,A)
%LORENZSTOP  Find local maxima of norm(y).
%  d/dt(norm(y)^2) = d/dt(y'*y) = 2*ydot'*y

A(1,3) = y(2);
A(3,1) = -y(2);
ydot = A*y;
gstop = ydot'*y;
isterm = 0;
direct = -1;
