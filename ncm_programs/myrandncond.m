% RANDNCOND  Condition of random matrices

nmax = 100;
n = 2:nmax;
p = 3/2;
kappalo = n.^p;
kappahi = 50*n.^p;

shg
clf reset
h = loglog(n,[kappalo; kappahi],'-',nmax,NaN,'.');
set(h(1:2),'color',[0 .5 0]);
set(gca,'xtick',[2:2:10 20:20:nmax])
kappamax = 1.e6;
axis([2 nmax 2 kappamax])
stop = uicontrol('pos',[20 10 40 25], ...
   'style','toggle','string','stop','value',0);

h = h(3);
set(h,'erasemode','none','color','blue')
while get(stop,'value') ~= 1
   n = ceil(rand*nmax);
   A = randn(n,n);
   kappa = condest(A,1);
   set(h,'xdata',n,'ydata',kappa)
   drawnow
end
