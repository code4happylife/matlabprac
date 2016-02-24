ms = 16;
format long e

disp('Norris')
load norris.dat
y = norris(:,1);
x = norris(:,2);
A = [ones(size(x)) x];
beta = A\y
u = (0:1000)';
v = beta(1)+beta(2)*u;
subplot(2,2,1)
plot(x,y,'.',u,v,'-','markersize',ms);
axis([-40 1040 -40 1040])
ylabel('y')
title('Norris')

subplot(2,2,3)
r = y - (beta(1)+beta(2)*x);
plot(x,r,'.','markersize',ms)
line([-40 1040],[0 0],'color','k')
axis([-40 1040 -2.5 2.5])
xlabel('x')
ylabel('r')

disp('Pontius')
load pontius.dat
y = pontius(:,1);
x = pontius(:,2);
A = [ones(size(x)) x x.^2];
beta = A\y
u = (1.5e5:.5e5:30e5)';
v = beta(1)+beta(2)*u+beta(3)*u.^2;

subplot(2,2,2)
plot(x,y,'.',u,v,'-','markersize',ms)
axis([0 3.2e6 0 2.5])
title('Pontius')

subplot(2,2,4)
r = y - (beta(1)+beta(2)*x+beta(3)*x.^2);
plot(x,r,'.','markersize',ms)
line([0 3.2e6],[0 0],'color','k')
axis([0 3.2e6 -5e-4 5e-4])
xlabel('x')
