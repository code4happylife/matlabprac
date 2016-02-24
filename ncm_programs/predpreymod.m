% PREDPREYMOD.  Modified Lotka-Volterra.

LV  = @(t,y) [2*y(1)-.01*y(1)*y(2); -y(2)+.01*y(1)*y(2)];
MLV = @(t,y) [2*(1-y(1)/400)*y(1)-.01*y(1)*y(2); -y(2)+.01*y(1)*y(2)];
[t,y] = ode23tx(LV,[0 50],[300 150]');
[T,Y] = ode23tx(MLV,[0 50],[300 150]');

subplot(2,2,1)
plot(t,y)
axis([0 50 0 500])
legend('rabbits','foxes')
title('LV')
xlabel('time')

subplot(2,2,2)
plot(T,Y)
axis([0 50 0 500])
legend('rabbits','foxes')
title('Lotka-Volterra')
xlabel('time')

subplot(2,2,3)
plot(y(:,2),y(:,1))
axis([0 480 0 360])
axis equal
title('LV')
xlabel('foxes')
ylabel('rabbits')

subplot(2,2,4)
plot(Y(:,2),Y(:,1))
axis([0 480 0 360])
axis equal
title('Modified Lotka-Volterra')
xlabel('foxes')
ylabel('rabbits')
