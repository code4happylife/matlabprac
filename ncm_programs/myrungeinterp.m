function rungeinterp(arg)
%RUNGEINTERP  Runge's polynomial interpolation example.
%   F(x) = 1/(1+25*x^2)
%   Polynomial interpolation at equally spaced points, -1 <= x <= 1.
%   Does interpolant converge as number of points is increased?

F = @(x) 1./(1+25*x.^2);

if nargin == 0

   % Initialize plot and uicontrols

   shg
   clf reset
   set(gcf,'doublebuffer','on','numbertitle','off', ...
       'name','Runge''s interpolation example')
   n = 1;
   u = -1.1:.01:1.1;
   z = F(u);
   h.plot = plot(u,z,'-', 0,1,'o', u,z,'-', [0 0],[1 1],'*');
   set(h.plot(1),'color',[.6 .6 .6]);
   set(h.plot(2),'color','blue');
   set(h.plot(3),'color',[0 2/3 0]);
   set(h.plot(4),'color','black');
   axis([-1.1 1.1 -0.1 1.1])
   title(char(F),'interpreter','none')

   h.minus = uicontrol('units','norm','pos',[.38 .01 .06 .05], ...
          'fontsize',12,'string','<','callback','myrungeinterp(''n--'')');
   h.n = uicontrol('units','norm','pos',[.46 .01 .12 .05], ...
          'fontsize',12,'userdata',n,'callback','myrungeinterp(''n=1'')');
   h.plus = uicontrol('units','norm','pos',[.60 .01 .06 .05], ...
          'fontsize',12,'string','>','callback','myrungeinterp(''n++'')');
   h.close = uicontrol('units','norm','pos',[.80 .01 .10 .05], ...
          'fontsize',12,'string','close','callback','close');
   h.cheby = uicontrol('units','norm','pos',[.05 .01 .25 .05], ...
          'fontsize',12,'string','Chebyshev points','style','toggle', ...
          'callback','myrungeinterp(''n'')');
   h.star = text(.3,-.01,'*');

   set(gcf,'userdata',h)
   arg = 'n=1';
end

% Update plot.

h = get(gcf,'userdata');

% Number of interpolation points.

n = get(h.n,'userdata');
switch arg
   case 'n--', n = n-2;
   case 'n++', n = n+2;
   case 'n=1', n = 1;
   case 'n', % Leave n unchanged.
end
set(h.n,'string',['n = ' num2str(n)],'userdata',n);
if n==1
   set(h.minus,'enable','off');
else
   set(h.minus,'enable','on');
end

if n == 1;
   x = 0;
else
   x = -1 + 2*(0:n-1)/(n-1);
   if get(h.cheby,'val') == 1
      x = sin(pi/2*x);
   end
end
y = F(x);
u = get(h.plot(1),'xdata');
v = polyinterp(x,y,u);
k = min(find(u>0 & abs(v-F(u))>.025));
k = max(find(x<u(k)));
ustar = (x(k)+x(min(n,k+1)))/2;
vstar = polyinterp(x,y,ustar);
set(h.plot(2),'xdata',x,'ydata',y);
set(h.plot(3),'xdata',u,'ydata',v);
set(h.plot(4),'xdata',[-ustar ustar],'ydata',[vstar vstar]);
set(h.star,'string',sprintf('x_* = %6.4f',ustar))
