function peakspdes
% PEAKSPDES   PDES using PEAKS

peaks_poisson
toggle('next')
peaks_heat_source
toggle('next')
peaks_heat_initial_val
toggle('next')
peaks_wave
toggle('close')
close(gcf)


% ------------------------

function peaks_poisson(n);
% 2-D Poisson equation, source = peaks(x,y), boundary values = 0.

if nargin < 1, n = 31; end
h = 6/(n+1);
[x,y] = meshgrid(-3:h:3);
f = peaks(x,y);
e = ones(n,1);
D = spdiags([e -2*e e],[-1 0 1],n,n);
I = eye(n,n);
A = (1/h^2)*(kron(D,I) + kron(I,D));
j = 2:n+1;
F = f(j,j);
U = A\F(:);
u = zeros(n+2,n+2);
u(j,j) = reshape(U,n,n);

% subplot(2,2,1)
% contourf(x,y,f,-8:1:8)
% caxis([-8 8])
% axis square
% set(gca,'ytick',get(gca,'xtick'))

% subplot(2,2,2)
a = min(min(u));
b = max(max(u));
c = a + (b-a)*(0:.1:10)';
contourf(x,y,u,c)
caxis([a b])
% axis square
colorbar
set(gca,'ytick',get(gca,'xtick'))


% ------------------------

function peaks_heat_source(n,sigma)
% 2-D heat equation
% initial values = 0, source and boundary values = -peaks(x,y)

if nargin < 1, n = 31; end
if nargin < 2, sigma = 1/4; end
set(gcf,'doublebuffer','on')
h = 6/(n+1);
[x,y] = meshgrid(-3:h:3);
u = zeros(size(x));
f = h^2*peaks(x,y);
t = 0;
delta = sigma*sqrt(h);
p = 2:n+1;
stop = uicontrol('units','norm','pos',[.02 .02 .10 .04], ...
   'string','stop','style','toggle');
while ~get(stop,'value')
   t = t + delta;
   w = u;
   u(p,p) = u(p,p) + ...
      sigma*(u(p+1,p) + u(p-1,p) + u(p,p+1) + u(p,p-1) - 4*u(p,p) - f(p,p));
   a = min(min(u));
   b = max(max(u));
   c = a + (b-a)*(0:.1:1);
   contourf(x,y,u,c)
   caxis([a b])
   colorbar
   title(sprintf('%10.4f',t))
   drawnow
end


% ------------------------

function peaks_heat_initial_val(n,sigma)
% 2-D heat equation
% initial values = peaks(x,y), source and boundary values = zero.

if nargin < 1, n = 31; end
if nargin < 2, sigma = 1/4; end
set(gcf,'doublebuffer','on')
h = 6/(n+1);
[x,y] = meshgrid(-3:h:3);
u = peaks(x,y);
t = 0;
delta = sigma*sqrt(h);
p = 2:n+1;
stop = uicontrol('units','norm','pos',[.02 .02 .10 .04], ...
   'string','stop','style','toggle');
while ~get(stop,'value')
   t = t + delta;
   w = u;
   u(p,p) = u(p,p) + ...
      sigma*(u(p+1,p) + u(p-1,p) + u(p,p+1) + u(p,p-1) - 4*u(p,p));
   a = min(min(u));
   b = max(max(u));
   c = a + (b-a)*(0:.1:1);
   contourf(x,y,u,c)
   caxis([a b])
   colorbar
   title(sprintf('%10.4f',t))
   drawnow
end


% ------------------------

function peaks_wave(n,sigma);
% 2-D wave equation
% initial values = peaks(x,y)
% initial velocity, source and boundary values = zero.

set(gcf,'doublebuffer','on')
if nargin < 1, n = 63; end
h = 6/(n+1);
[x,y] = meshgrid(-3:h:3);
u = peaks(x,y);
contourf(x,y,u,-8:1:8)
caxis([-8 8])
colorbar
t = 0;
title(sprintf('%8.4f',t))

titl = title(num2str(t));
if nargin < 2, sigma = 1/2; end
delta = sqrt(sigma)*h;
v = u;
w = u;
p = 2:n+1;
stop = uicontrol('units','norm','pos',[.02 .02 .10 .04], ...
   'string','stop','style','toggle');
while ~get(stop,'value')
   t = t + delta;
   w(p,p) = 2*u(p,p) - v(p,p) + ...
      sigma*(u(p+1,p) + u(p-1,p) + u(p,p+1) + u(p,p-1) - 4*u(p,p));
   v = u;
   u = w;
   contourf(x,y,u,-8:1:8)
   caxis([-8 8])
   colorbar
   title(sprintf('%10.4f',t))
   drawnow
end


% ------------------------

function toggle(str)
uic = uicontrol('units','norm','pos',[.02 .02 .10 .04], ...
   'string',str,'style','toggle');
while get(uic,'value') == 0
   drawnow
end
