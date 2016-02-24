% INTERPS  Four interpolants.
shg
clf reset
x = 1:6;
y = [16 18 21 17 15 12];
u = .75:.05:6.25;
c = get(gca,'colororder');

subplot(2,2,1)
h = plot(x,y,'o',x,y,'-');
set(h(1),'color','black')
set(h(2),'color',c(2,:))
axis([0.5 6.5 10 22])
title('Piecewise linear interpolation')
set(gca,'xtick',[],'ytick',[])

subplot(2,2,2)
p = polyinterp(x,y,u);
h = plot(x,y,'o',u,p,'-');
set(h(1),'color','black')
set(h(2),'color',c(3,:))
axis([0.5 6.5 10 22])
title('Full degree polynomial interpolation')
set(gca,'xtick',[],'ytick',[])

subplot(2,2,3)
q = pchip(x,y,u);
h = plot(x,y,'o',u,q,'-');
set(h(1),'color','black')
set(h(2),'color',c(5,:))
axis([0.5 6.5 10 22])
title('Shape preserving Hermite interpolation')
set(gca,'xtick',[],'ytick',[])

subplot(2,2,4)
s = splinetx(x,y,u);
h = plot(x,y,'o',u,s,'-');
set(h(1),'color','black')
set(h(2),'color',c(4,:))
axis([0.5 6.5 10 22])
title('Spline interpolation')
set(gca,'xtick',[],'ytick',[])
