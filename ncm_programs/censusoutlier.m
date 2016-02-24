function censusgui(arg)
%CENSUSGUI Try to predict the US population in the year 2010.
% This example is older than MATLAB.  It started as an exercise in
% "Computer Methods for Mathematical Computations", by Forsythe,
% Malcolm and Moler, published by Prentice-Hall in 1977.
% The data set has been updated every ten years since then.
% Today, MATLAB makes it easier to vary the parameters and see the
% results, but the underlying mathematical principles are unchanged:
%
%    Using polynomials of even modest degree to predict
%    the future by extrapolating data is a risky business.
%
% The data is from the decennial census of the United States for the
% years 1900 to 2000.  The task is to extrapolate beyond 2000.
% In addition to polynomials of various degrees, you can choose
% interpolation by a cubic spline, interpolation by a shape-preserving
% Hermite cubic, and a least squares fit by an exponential.
% Error estimates attempt to account for errors in the data,
% but not in the extrapolation model.

% Census data for 1900 to 2000.
% The population on April 1, 2000 was 281,421,906, according to:
% http://www.census.gov/main/www/cen2000.html

p = [ 75.995  91.972 105.711 123.203 131.669 150.697 ...
     179.323 203.212 226.505 249.633 281.422]';

% Create outlier
p(6) = 50.697;

t = (1900:10:2000)';   % Census years
x = (1890:1:2019)';    % Evaluation years
w = 2010;              % Extrapolation target
guess = 320;           % Wild guess
dmax = length(t)-1;    % Maximum polynomial degree

if nargin == 0

   % Initialize plot and uicontrols

   shg
   clf reset
   set(gcf,'doublebuffer','on','name','Census gui', ...
       'menu','none','numbertitle','off')
   h.plot = plot(t,p,'bo', x,0*x,'k-', w,0,'.', [x;NaN;x],[x;NaN;x],':');
   darkgreen = [0 2/3 0];
   marksize = get(0,'defaultlinemarkersize');
   set(h.plot(3:4),'color',darkgreen,'markersize',4*marksize-6);
   axis([min(x) max(x) 0 400]);
   title('Predict U.S. Population in 2010');
   ylabel('Millions');

   h.text = text(w-4,guess,'?','color',darkgreen,'fontweight','bold');
   h.model = uicontrol('units','norm','pos',[.20 .80 .20 .05], ...
           'style','popup','background','white','string', ...
           {'census data','polynomial','pchip','spline','exponential'}, ...
           'callback','censusoutlier(''cb'')');
   h.deg = uicontrol('units','norm','pos',[.20 .75 .20 .05], ...
           'style','text','background','white', ...
           'horiz','left','string','?');
   h.slide = uicontrol('units','norm');
   set(h.slide,'pos',[.20 .72 .20 .04])
   set(h.slide,'style','slider','background','white')
   set(h.slide,'min',0,'max',dmax)
   set(h.slide,'sliderstep',[1 2]/dmax)
   set(h.slide,'value',3)
   set(h.slide,'callback','censusoutlier(''cb'')');
   h.err = uicontrol('units','norm','pos',[.20 .65 .20 .05], ...
           'style','check','background','white', ...
           'string','error estimates', ...
           'callback','censusoutlier(''cb'')');
   set(gcf,'userdata',h);
   uicontrol('style','push','units','normal','pos',[.02 .01 .10 .06], ...
      'string','close','callback','close(gcf)')

end

% Update plot with new model
% y = interpolated values
% z = extrapolated value

h = get(gcf,'userdata');
models = get(h.model,'string');
model = models{get(h.model,'value')};
d = round(get(h.slide,'value'));  % Polynomial degree

switch model
   case 'census data'
      y = NaN*x;
      z = guess;
   case 'polynomial'
      s = (t-1950)/50;   c = polyfit(s,p,d);
      s = (x-1950)/50;   y = polyval(c,s);
      s = (w-1950)/50;   z = polyval(c,s);
   case 'pchip'
      y = pchip(t,p,x);
      z = pchip(t,p,w);
   case 'spline'
      y = spline(t,p,x);
      z = spline(t,p,w);
   case 'exponential'
      c = polyfit(log(t),log(p),1);
      y = exp(polyval(c,log(x)));
      z = exp(polyval(c,log(w)));
end
set(h.plot(2),'ydata',y);
set(h.plot(3),'ydata',z);
set(h.text,'pos',[w-18,min(max(z,20),380)],'string',sprintf('%8.3f',z))

switch model
   case 'census data'
      set(h.err,'vis','off','value',0);
      set([h.deg; h.slide],'vis','off');
      set(h.text,'pos',[w-4,z],'string','?')
   case 'polynomial'
      set(h.err,'vis','on','pos',[.20 .65 .20 .05]);
      set([h.deg; h.slide],'vis','on');
      set(h.deg,'string',['degree = ' num2str(d)]);
   otherwise
      set(h.err,'vis','on','pos',[.20 .75 .20 .05]);
      set([h.deg; h.slide],'vis','off');
end

if get(h.err,'value');
   errest = errorestimates(model,t,p,x,y,d);
   set(h.plot(4),'vis','on','ydata',errest);
else
   set(h.plot(4),'vis','off');
end


% ------------------------------------------------

function errest = errorestimates(model,t,p,x,y,d)
% Provide error estimates for censusgui

switch model
   case 'polynomial'
      if d > 0
         V(:,d+1) = ones(size(t));
         s = (t-1950)/50;
         for j = d:-1:1
            V(:,j) = s.*V(:,j+1);
         end
         [Q,R] = qr(V);
         R = R(1:d+1,:);
         RI = inv(R);
         E = zeros(length(x),d+1);
         s = (x-1950)/50;
         for j = 1:d+1
            E(:,j) = polyval(RI(:,j),s);
         end
         sig = 10;   % Rough estimate
         e = sig*sqrt(1+diag(E*E'));
         errest = [y-e; NaN; y+e];
      else
         errest = [y-NaN; NaN; y+NaN];
      end
   case {'pchip','spline'}
      n = length(t);
      I = eye(n,n);
      E = zeros(length(x),n);
      for j = 1:n
         if isequal(model,'pchip')
            E(:,j) = pchip(t,I(:,j),x);
         else
            E(:,j) = spline(t,I(:,j),x);
         end
      end
      sig = 10;  % Rough estimate
      e = sig*sqrt(1+diag(E*E'));
      errest = [y-e; NaN; y+e];
   case 'exponential'
      V = [ones(size(t)) log(t)];
      [Q,R] = qr(V);
      c = R\(Q'*log(p));
      r = log(p) - V*c;
      E = [ones(size(x)) log(x)]/R(1:2,1:2);
      sig = norm(r);
      e = sig*sqrt(1+diag(E*E'));
      errest = [y.*exp(-e); NaN; y.*exp(e)];
end
