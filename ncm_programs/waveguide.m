% WAVEGUIDE  The first eigenfunction of the H-shaped ridge waveguide.
shg
clf
set(gcf,'color','w')

h = 1/32;
[x,y] = meshgrid(-1:h:1);
xv = [0 1 1 2 2 1 1 -1 -1 -2 -2 -1 -1 0]/2; 
yv = [1 1 2 2 -2 -2 -1 -1 -2 -2 2 2 1 1]/2;
[in,on] = inregion(x,y,xv,yv);
p = find(in-on);
G = zeros(size(x));
G(p) = 1:length(p);
A = delsq(G);
n = size(A,1);
opts.disp = 0;
[u,lam] = eigs(A,1,0,opts);
m = find(abs(u) == max(abs(u)));
U = zeros(size(x));
p = find(G>0);
U(p) = u/u(m);

contour(x,y,U,[.05:.1:1])
H = hot(18);
colormap(H(2:11,:))
axis square
axis off

x = [0 0 1];
y = [1/2 0 0];
line(xv,yv,'color','k','linestyle','-','linewidth',2)
line(x,y,'color','k','linestyle','--','linewidth',2)
