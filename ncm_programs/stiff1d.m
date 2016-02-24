% STIFF1D.  Mildly stiff one dimensional equation.

F = @(t,y) -1000*(y-sin(t))+cos(t)
y = dsolve('Dy = -1000*(y-sin(t))+cos(t)','y(0)=1')

[t,y] = ode23tx(F,[0 1],1);
nsteps23tx = length(t)-1

[T,Y] = ode23s(F,[0 1],1);
nsteps23s = length(T)-1

subplot(2,1,1)
h = plot(T,Y,'o',t,y,'-');
set(h(1),'markersize',8)

subplot(2,2,3)
h = plot(T,Y,'o',t,y,'.');
axis([0 .02 .005 .02])
set(h(1),'markersize',8)
set(h(2),'markersize',16)

subplot(2,2,4)
h = plot(T,Y,'o',t,y,'.');
axis([.4 .55 .38 .53])
set(h(1),'markersize',8)
set(h(2),'markersize',16)
