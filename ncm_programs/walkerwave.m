function walker
% WALKER  Human gait.
% This model, developed by Nikolaus Troje, is a five-term Fourier series
% with vector-valued coefficients that are the principal components for
% data obtained in motion capture experiments involving subjects wearing
% reflective markers walking on a treadmill.  The components, which are
% also known as "postures" or "eigenwalkers", correspond to the static
% position, forward motion, sideways sway, and two hopping/bouncing
% movements that differ in the phase relationship between the upper and
% lower portions of the body.  The postures are also classified by gender.
% Sliders allow you to vary the amount that each component contributes to
% the overall motion.  A slider setting greater than 1.0 overemphasizes
% the characteristic.  Can you see whether positive values of the gender
% coefficient correspond to male or female subjects?
%
% References:
%    http://www.bml.psy.ruhr-uni-bochum.de/Demos
%    http://www.biomotionlab.de/Text/WDP2002_Troje.pdf
%    http://journalofvision.org/2/5/2

clf
shg
set(gcf,'doublebuf','on','color','w','name','Walker','numbertitle','off')
set(gca,'pos',get(gca,'pos')+[0 .07 0 0])

% The body is represented by 15 points in three space, i.e. a vector of
% length 45.  The data consists of F, five vectors describing the average
% female and M, five vectors describing the average male.  Three linked
% segments, indexed by L, are the head/torso, the arms, and the legs.

% Initial view

load walkers
xyz = reshape((F(:,1)+M(:,1))/2,15,3);
L = {[1 5 12],[2 3 4 5 6 7 8],[9 10 11 12 13 14 15]};
for k = 1:3
   p(k) = line(xyz(L{k},1),xyz(L{k},2),xyz(L{k},3),'marker','o', ...
      'markersize',10,'linestyle','-','erasemode','background');
end
axis([-750 750 -750 750 0 1500])
set(gca,'xtick',[],'ytick',[],'ztick',[])
view(160,10)
cameratoolbar
uicontrol('style','text','units','norm','pos',[.00 .92 .25 .06], ...
   'back','white','fontangle','italic','string', ...
   {'Change the view','with the mouse'})

% Sliders and controls

labels = {'speed','gender','stride','sway','hop','bounce'};
for j = 1:6
   switch j
      case 1, smin = 0; start = 0.5; smax = 2;
      case 2, smin = -3; start = 0; smax = 3;
      otherwise, smin = -2; start = 1; smax = 2;
   end
   txt = uicontrol('style','text','string',sprintf('%4.2f',start), ...
      'back','w','units','norm','pos',[.16*j-.10 .10 .08 .03]);
   sliders(j) = uicontrol('style','slider','units','norm','back','w', ...
      'pos',[.16*j-.13 .06 .14 .03],'min',smin,'max',smax,'val',start, ...
      'sliderstep',[1/(4*smax),1/(10*smax)],'userdata',txt,'callback',...
      'set(get(gco,''userd''),''str'',sprintf(''%4.2f'',get(gco,''val'')))');
   uicontrol('style','text','string',labels{j},'back','w', ...
      'units','norm','pos',[.16*j-.10 .02 .08 .03])
end
txt = uicontrol('style','text','string','0.00', ...
   'back','w','units','norm','pos',[.06 .26 .08 .03]);
sliders(7) = uicontrol('style','slider','units','norm','back','w', ...
   'pos',[.03 .22 .14 .03],'min',0,'max',1,'val',0, ...
   'sliderstep',[1/10,1/20],'userdata',txt,'callback',...
   'set(get(gco,''userd''),''str'',sprintf(''%4.2f'',get(gco,''val'')))');
uicontrol('style','text','string','wave','back','w', ...
   'units','norm','pos',[.06 .18 .08 .03])
stop = uicontrol('style','toggle','units','norm','pos',[.92 .94 .08 .06], ...
   'back','w','string','stop');

% Start walkin'...

period = 151.5751;
omega = 2*pi/period;
t = 0;
while get(stop,'value') == 0
   s = cell2mat(get(sliders,'value'));
   t = t + s(1);
   c = [sin(omega*t); cos(omega*t); sin(2*omega*t); cos(2*omega*t)];
   w = [1; s(3:6).*c];
   f = reshape(F*w,15,3);
   m = reshape(M*w,15,3);
   xyz = (f+m)/2 + s(2)*(f-m)/2;
   theta1 = s(7)*pi/2*(1+.25*cos(2*omega*t));
   theta2 = s(7)*pi/4;
   wave = [ 0  200*cos(theta1) 500*sin(theta1)
            0  100*cos(theta2) 250*sin(theta2)];
   xyz(2:3,:) = xyz(2:3,:) + wave;
   for k = 1:3
      set(p(k),'xdata',xyz(L{k},1),'ydata',xyz(L{k},2),'zdata',xyz(L{k},3));
   end
   drawnow
end;
set(stop,'value',0,'string','close','fontweight','bold','callback','close(gcf)')
