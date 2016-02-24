function ppp(r0,f0)
% Periodic preditor prey

y0 = [r0; f0];
alpha = .01;
opts = odeset('reltol',1.e-6,'abstol',1.e-6, ...
   'refine',4,'events',@g);
[t,y,te,ye] = ode45(@f,[0 25],y0,opts,alpha,y0);
te
ye

shg
subplot(2,2,1)
plot(t,y)
xlabel(['tp = ' num2str(t(end))])

subplot(2,2,2)
plot(r0,f0,'o',y(:,1),y(:,2),'-')
xlabel('rabbits')
ylabel('foxes')

% ------------------------------------------

function ydot = f(t,y,alpha,y0)
ydot = [2*y(1)-alpha*y(1)*y(2)
         -y(2)+alpha*y(1)*y(2)];


% ------------------------------------------

function [gstop,isterm,dir] = g(t,y,alpha,y0)
ydot = f(t,y,alpha,y0);
d = (y-y0)./y;
gstop = d'*ydot;
isterm = 1;
dir = 1;
