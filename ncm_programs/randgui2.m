function p = randgui2(randfun);
% RANDGUI2. Monte Carlo computation of four times the ratio of
% the area of a circle to the area of a surrounding square.

if nargin < 1, randfun = @rand; end
figure(gcf)
clf reset
set(gcf,'doublebuffer','on')
z = exp(2*pi*i*(0:255)/256);
plot(real(z),imag(z),'r-','linewidth',2)
axis square
set(gca,'xticklabel',[],'yticklabel',[])
stop = uicontrol('style','toggle','string','stop','value',0);
t = text(1.1,0.0,num2str(pi),'fontsize',16);
hold on
k = 0;
n = 0;
while get(stop,'value') == 0
   x = 2*randfun()-1;
   y = 2*randfun()-1;
   plot(x,y,'b.','markersize',16)
   drawnow
   if x^2 + y^2 < 1
      k = k+1;
   end
   n = n+1;
   p = 4*k/n;
   set(t,'string',num2str(p))
end
set(stop,'string','close','value',0,'callback','close(gcf)')
