function humpspdes(n)
if nargin < 1, n = 31; end

humps_poisson(2*n+1)
toggle('next')
humps_heat_source(n)
toggle('next')
humps_heat_initial_val(n)
toggle('next')
humps_wave(4*n+3)
toggle('close')
close(gcf)

% ------------------------------

function humps_poisson(n);
% 1-D Poisson equation, source = -humps(x), boundary values = 0.
if nargin < 1, n = 63; end
shg
clf
h = 1/(n+1);
x = (0:h:1)';
f = -humps(x);
e = ones(n,1);
u = zeros(n+2,1);
A = spdiags([e -2*e e],[-1 0 1],n,n)/h^2;
j = 2:n+1;
u(j) = A\f(j);
subplot(2,2,1)
plot(x,f)
subplot(2,2,2)
plot(x,u);
axis([0 1 -1 5])
subplot(2,2,3)
residual = diff(u,2)/h^2 + humps(x(j));
plot(x(j),residual)


% ------------------------------

function humps_heat_source(n,sigma);
% 1-D heat equation, source = humps(x), initial and boundary values = 0.
if nargin < 1, n = 31; end
if nargin < 2, sigma = 1/2; end
shg
clf
set(gcf,'doublebuffer','on')
h = 1/(n+1);
x = (0:h:1)';
u = zeros(n+2,1);
f = h^2*humps(x);
p = plot(x,u);
axis([0 1 -1 5])
t = 0;
titl = title(num2str(t));
delta = h^2*sigma;
j = 2:n+1;
stop = uicontrol('units','norm','pos',[.02 .02 .10 .04], ...
   'string','stop','style','toggle');
while get(stop,'value') == 0 & t < 1
   t = t + delta;
   u(j) = u(j) + sigma*(u(j+1) - 2*u(j) + u(j-1) + f(j));
   set(p,'ydata',u)
   set(titl,'string',sprintf('%8.4f',t))
   drawnow
end


% ------------------------------

function humps_heat_initial_val(n,sigma);
% 1-D heat equation, source = 0, initial and boundary values = humps(x).
if nargin < 1, n = 31; end
if nargin < 2, sigma = 1/2; end
shg
clf
set(gcf,'doublebuffer','on')
h = 1/(n+1);
x = (0:h:1)';
u = humps(x);
p = plot(x,u);
axis([0 1 -1.1*min(u) 1.1*max(u)])
t = 0;
titl = title(num2str(t));
delta = h^2*sigma;
j = 2:n+1;
stop = uicontrol('units','norm','pos',[.02 .02 .10 .04], ...
   'string','stop','style','toggle');
while get(stop,'value') == 0 & t < 1
   t = t + delta;
   u(j) = u(j) + sigma*(u(j+1) - 2*u(j) + u(j-1));
   set(p,'ydata',u)
   set(titl,'string',sprintf('%8.4f',t))
   drawnow
end


% ------------------------------

function humps_wave(n,sigma);
% 1-D wave equation, source = 0, initial and boundary values = humps(x).
if nargin < 1, n = 255; end
if nargin < 2, sigma = 1; end
shg
clf
set(gcf,'doublebuffer','on')
h = 1/(n+1);
x = (0:h:1)';
u = humps(x);
v = u;
w = u;
p = plot(x,u);
s = 1.1*max(abs(u));
axis([0 1 -s s])
t = 0;
titl = title(num2str(t));
delta = sqrt(sigma)*h;
j = 2:n+1;
stop = uicontrol('units','norm','pos',[.02 .02 .10 .04], ...
   'string','stop','style','toggle');
while get(stop,'value') == 0 & t < 4
   t = t + delta;
   w(j) = 2*u(j) - v(j) + sigma*(u(j+1) - 2*u(j) + u(j-1));
   v = u;
   u = w;
   set(p,'ydata',u)
   set(titl,'string',sprintf('%8.4f',t))
   drawnow
end


% ------------------------

function toggle(str)
uic = uicontrol('units','norm','pos',[.02 .02 .10 .04], ...
   'string',str,'style','toggle');
while get(uic,'value') == 0
   drawnow
end
