function sierpinski
%SIERPINSKI  Sierpinski's triangle.
%   This version runs forever, or until stop is toggled.
%   See also FINITESIERPINSKI.

shg
clf reset
set(gcf,'color','white','menubar','none', ...
   'numbertitle','off','name','Sierpinski''s triangle')
x = [0; 0];
h = plot(x(1),x(2),'.');
darkgreen = [0 2/3 0];
set(h,'markersize',1,'color',darkgreen,'erasemode','none');
axis([0 1 0 1])
axis square
axis off
stop = uicontrol('style','toggle','string','stop', ...
   'background','white');
drawnow

% Define three affine transformations, x = A*x + bk.

A = [1/2 0; 0 1/2];
b1 = [0  0]';
b2 = [1/2  0]';
b3 = [1/4  sqrt(3)/4]';

% Run forever

cnt = 1;
tic
while ~get(stop,'value')
   r = ceil(3*rand);
   switch r
      case 1, x = A*x + b1;
      case 2, x = A*x + b2;
      case 3, x = A*x + b3;
   end
   set(h,'xdata',x(1),'ydata',x(2));
   drawnow
   cnt = cnt + 1;
end
t = toc;
s = sprintf('%8.0f points in %6.3f seconds',cnt,t);
text(-1.5,-0.5,s,'fontweight','bold');
set(stop,'style','pushbutton','string','close','callback','close(gcf)')
