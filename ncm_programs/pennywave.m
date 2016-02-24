function pennywave(delta)
% PENNYWAVE
% Initial value of the height is obtained from measurements
% made at the National Institute of Science and Technology
% of the depth of a mold for a U. S. one cent coin.
% What is the limiting value of the height as t -> inf ?
% pennywave(delta) takes time steps of size delta.
% For what values of delta is the computation stable?

shg
clf
set(gcf,'doublebuffer','on')
% set(gcf,'renderer','zbuffer')
if nargin < 1, delta = .5; end
load penny
U = flipud(P);
V = U;
W = U;
surfu = surf(U);
daspect([1,1,128])
colormap(copper)
shading interp
material metal
lighting gouraud
view(2)
axis tight
axis off
set(gca,'zlimmode','auto','climmode','manual');
light('pos',[1,2,3000],'style','inf');
titl = title('0');
pause(2)
[p,q] = size(P);
n = [2:p p];
e = [2:q q];
s = [1 1:p-1];
w = [1 1:q-1];
p = 2:size(P,1)-1;
q = 2:size(P,2)-1;
t = 0;
h = 1;
sigma = (delta/h)^2;
stop = uicontrol('units','norm','pos',[.02 .02 .10 .05], ...
   'string','stop','style','toggle');
%while ~get(stop,'value')
while t < 15
   W = 2*U - V + sigma*(U(n,:)+U(s,:)+U(:,e)+U(:,w)-4*U);
   V = U;
   U = W;
   t = t + delta;
   set(surfu,'zdata',U,'cdata',U)
   set(titl,'string',sprintf('%10.3f',t))
   drawnow
end
set(stop,'str','close','val',0,'callback','close(gcf)')
