clf reset
shg

x = [1.02 .95 .87 .77 .67 .56 .44 .30 .16 .01]';
y = [ .39 .32 .27 .22 .18 .15 .13 .12 .13 .15]';
n = length(x);
plot(x,y,'b.','markersize',16)
hold on

A = [x.^2 x.*y y.^2 x y];
f = ones(n,1);
c = A\f
z = c(1)*x.^2 + c(2)*x.*y + c(3)*y.^2 + c(4)*x + c(5)*y - f;

[X,Y] = meshgrid(-1.25:.05:1.25);
Z = c(1)*X.^2 + c(2)*X.*Y + c(3)*Y.^2 + c(4)*X + c(5)*Y - 1;
contour(X,Y,Z,[0 0])
colormap(1-white)

x = x + .005*(2*rand(n,1)-1);
y = y + .005*(2*rand(n,1)-1);
plot(x,y,'o','markersize',12,'color',[0 2/3 0])
A = [x.^2 x.*y y.^2 x y];
c = A\f
Z = c(1)*X.^2 + c(2)*X.*Y + c(3)*Y.^2 + c(4)*X + c(5)*Y - 1;
contour(X,Y,Z,[0 0])
colormap(1-white)
axis equal
hold off
