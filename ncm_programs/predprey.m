function  e = predprey(r0,f0,tf)
% Lotka-Volterra

y0 = [r0; f0];
opts = odeset('reltol',1.e-5,'stats','on');
alpha = .01;
[t,y] = ode23(@f,[0 tf],y0,opts,alpha);

close all
figure(1)
plot(t,y)
xlabel('t')
legend('rabbits','foxes')

figure(2)
plot(y(1,1),y(1,2),'bo',y(end,1),y(end,2),'ro',y(:,1),y(:,2),'-')
xlabel('rabbits')
ylabel('foxes')

e = y(end,:)-y0';

% ------------------------------------------

function ydot = f(t,y,alpha)
ydot = [2*y(1)-alpha*y(1)*y(2)
         -y(2)+alpha*y(1)*y(2)];
