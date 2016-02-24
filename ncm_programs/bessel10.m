% Bessel10

x = 0:pi/50:10*pi;
bj = besselj(0,x);
bjz = zeros(1,10);
by = bessely(0,x);
byz = zeros(1,10);
e = zeros(1,10);
J0 = @(x) besselj(0,x);
Y0 = @(x) bessely(0,x);
JY0 = @(x) besselj(0,x)-bessely(0,x);
for k = 1:10
   bjz(k) = fzerotx(J0,[k-.5 k+.5]*pi);
   byz(k) = fzerotx(Y0,[k-1 k]*pi);
   e(k) = fzerotx(JY0,[k-1 k]*pi);
end
f = besselj(0,e);
set(gcf,'defaultaxescolororder', ...
   [0 .5 0; 0 0 1; 0 0 0; 0 .5 .5; 0 .5 1; 0 0 0])
plot(x,bj,'-', x,by,'-', [0 10*pi],[0 0],'-', ...
    bjz,0*bjz,'o', byz,0*byz,'o', e,f,'*');
axis([0 10*pi -.5 1])

set(get(gca,'children'),'markersize',8)

