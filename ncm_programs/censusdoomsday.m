function censusdoomsday(callbackarg)

p = [ 75.995  91.972 105.711 123.203 131.669 150.697 ...
     179.323 203.212 226.505 249.633 281.422]';

t = (1900:10:2000)';   % Census years
x = (1890:1:2019)';    % Evaluation years
dmax = length(t)-1;    % Maximum polynomial degree

if nargin == 0

   % Initialize plot and uicontrols

   shg
   clf reset
   set(gcf,'doublebuffer','on','name','Census gui', ...
       'menu','none','numbertitle','off')
   h.plot = plot(t,p,'bo', x,0*x,'k-');
   axis([min(x) max(x) 0 400])
   title('Polynomial Extrapolation Predicts Disaster')
   ylabel('Millions')

   h.deg = uicontrol('units','norm','pos',[.26 .75 .13 .04], ...
           'tag','degree','style','text','background','white', ...
           'userdata',3,'string','degree = 8');
   h.ls = uicontrol('units','norm','pos',[.20 .75 .05 .04], ...
           'style','push','string','<','fontweight','bold', ...
           'callback','censusdoomsday(''<'')');
   h.gt = uicontrol('units','norm','pos',[.40 .75 .05 .04], ...
           'style','push','string','>','fontweight','bold', ...
           'callback','censusdoomsday(''>'')');
   set(gcf,'userdata',h);
   uicontrol('style','push','units','normal','pos',[.02 .01 .10 .06], ...
      'string','close','callback','close(gcf)')
   callbackarg = [];

else

   h = get(gcf,'userdata');
 
end

% Polynomial degree

d = get(h.deg,'userdata');
if isequal(callbackarg,'<'), d = d - 1; end
if isequal(callbackarg,'>'), d = d + 1; end
set(h.deg,'userdata',d,'string',sprintf('degree = %d',d))
set([h.ls; h.gt],'enable','on')
if d == 0, set(h.ls,'enable','off'), end
if d == dmax, set(h.gt,'enable','off'), end

% Case 'polynomial'

s = (t-1950)/50;
c = polyfit(s,p,d);
y = polyval(c,(x-1950)/50);
set(h.plot(2),'ydata',y);

% Find when population becomes zero.

if any(y(x > 2000) < 0)
   r = roots(c);
   r = r(imag(r)==0 & real(r)>0);
   z = 1950+50*r;
   date = datestr(datenum(z,4,1),22);
   text(z-4,-15,date);
   text(z-1,0,'X','color','r','fontweight','bold')
end
