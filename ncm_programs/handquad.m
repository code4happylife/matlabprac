function handquad
% HANDQUAD

% Load data.

load myhand.dat
x = myhand(:,1);
y = myhand(:,2);
ax = [-.05 1.05 -.05 1.05];

% Area of polygon

shg
clf
subplot(3,3,1)
Q = (x'*y([2:end 1]) - x([2:end 1])'*y)/2
p = plot(x,y,'bo',x,y,'k-');
set(p(1),'markersize',3)
set(p(2:end),'linewidth',1)
axis(ax)
set(gca,'xtick',[],'ytick',[])
title(sprintf('Q = %6.4f',Q))
drawnow

% Grid

subplot(3,3,2)
h = .05;
[u,v] = meshgrid(ax(1):h:ax(2),ax(3):h:ax(4));
k = inpolygon(u+h/2,v+h/2,x,y);
z = NaN*ones(size(u));
z(find(k)) = 0;
pcolor(u,v,z);
axis(ax)
hold on
p = plot(x,y,'bo',x,y,'k-');
set(p(1),'markersize',3)
hold off
Q = nnz(k)*h^2
set(gca,'xtick',[],'ytick',[])
title(sprintf('Q = %6.4f',Q))

% Numerical quadrature in 2-D

subplot(3,3,3)
p = plot(x,y,'bo',x,y,'k-');
set(p(1),'markersize',3)
axis(ax)
hold on
Q = dblquad(@chi,ax(1),ax(2),ax(3),ax(4),.001,[],x,y,p)
set(gca,'xtick',[],'ytick',[])
title(sprintf('Q = %6.4f',Q))
hold off

% print -depsc2 handquad.eps

% ------------------------

function k = chi(u,v,x,y,p)
if all(size(u) == 1), u = u(ones(size(v))); end
if all(size(v) == 1), v = v(ones(size(u))); end
k = inpolygon(u,v,x,y);
if nargin == 5
   plot(u(find(k)),v(find(k)),'.','markersize',6)
   drawnow
end
