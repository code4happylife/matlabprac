% CAPACITY  Electrostatic capacity of three regions

shg
clf
C = bone(128);
C = C(:,[2 3 1]);
colormap(C);

% Capacity of a square

m = 200;
h = 1/m;
S = numgrid('S',m+1);
A = 1/h^2 * delsq(S);
b = ones(nnz(S),1);
u = A\b;
U = reshape(u,m-1,m-1);
subplot(2,2,1)
contourf(U,10)
axis square
axis off
c = h^2*sum(sum(U))
title(sprintf('%7.5f',c))
drawnow

% Capacity of L-shaped domain

m = 200;
h = 1/m;
L = numgrid('L',2*m+1);
A = 1/h^2 * delsq(L);
b = ones(nnz(L),1);
u = A\b;
U = zeros(2*m+1,2*m+1);
U(find(L)) = u;
U(U==0) = NaN;
subplot(2,2,2)
contourf(U,10)
axis square
axis off
c = h^2*sum(u)
title(sprintf('%7.5f',c))
drawnow

% Capacity of my hand.

load myhand.dat
xv = myhand(:,1);
yv = myhand(:,2);
xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;
m = 200;
h = 1/m;
[x,y] = meshgrid(xmin:h:xmax, ymin:h:ymax);
[in,on] = inregion(x,y,xv,yv);
p = find(in-on == 1);
q = find(in-on == 0);
G = zeros(size(x));
G(p) = 1:length(p);
A = 1/h^2 * delsq(G);
b = ones(nnz(G),1);
u = A\b;
U = zeros(size(x));
U(p) = u;
U(q) = NaN;
axes('pos',[.354 .110 .327 .344])
contourf(U,10)
axis square
axis off
c = h^2*sum(u)
title(sprintf('%7.5f',c))
drawnow
