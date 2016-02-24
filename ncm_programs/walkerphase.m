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

labels = {'speed','stride','phase1','hop','phase2','gender'};
for j = 1:6
   switch j
      case 1, smin = 0; start = 0.5; smax = 2;
      case 6, smin = -3; start = 0; smax = 3;
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
stop = uicontrol('style','toggle','units','norm','pos',[.92 .94 .08 .06], ...
   'back','w','string','stop');

% Start walkin'...

F = reshape(F,15,3,5);
M = reshape(M,15,3,5);

fmean = F(:,:,1);
famp1 = sqrt(F(:,:,2).^2+F(:,:,3).^2);
fphi1 = atan2(F(:,:,3),F(:,:,2));
famp2 = sqrt(F(:,:,4).^2+F(:,:,5).^2);
fphi2 = atan2(F(:,:,5),F(:,:,4));

mmean = M(:,:,1);
mamp1 = sqrt(M(:,:,2).^2+M(:,:,3).^2);
mphi1 = atan2(M(:,:,3),M(:,:,2));
mamp2 = sqrt(M(:,:,4).^2+M(:,:,5).^2);
mphi2 = atan2(M(:,:,5),M(:,:,4));

period = 151.5751;
omega = 2*pi/period;
t = 0;
while get(stop,'value') == 0
   s = cell2mat(get(sliders,'value'));
   t = t + s(1);
   f = fmean + s(2)*famp1.*sin(omega*t+s(3)*fphi1) ...
             + s(4)*famp2.*sin(2*omega*t+s(5)*fphi2);
   m = mmean + s(2)*mamp1.*sin(omega*t+s(3)*mphi1) ...
             + s(4)*mamp2.*sin(2*omega*t+s(5)*mphi2);
   xyz = (f+m)/2 + s(6)*(f-m)/2;
   for k = 1:3
      set(p(k),'xdata',xyz(L{k},1),'ydata',xyz(L{k},2),'zdata',xyz(L{k},3));
   end
   drawnow
end;
set(stop,'value',0,'string','close','fontweight','bold','callback','close(gcf)')
