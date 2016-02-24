function biorhythm(birthday)
% When do all three biorhythms return to zero simultaneously?
if nargin < 1
   birthday = input('When is your birthday? ','s');
end

delta = lcm(lcm(23,28),33)/2;
y = fix(delta/365.25);
m = fix((delta - 365.25*y)/30.5);
d = delta - 365.25*y - 30.5*m;
disp(['delta = ' num2str(delta) ' days' ...
      ' = ' num2str(y) ' years + ' num2str(m) ' month + ' num2str(d) ' days'])

subplot(2,1,1)
t0 = datenum(birthday);
t1 = t0+delta;
t = (t1-28):1:(t1+28);
y = 100*[sin(2*pi*(t-t0)/23); sin(2*pi*(t-t0)/28); sin(2*pi*(t-t0)/33)];
plot(t,y)
line([t1 t1],[-100 100],'color','k')
set(gca,'xtick',(t1-28):7:(t1+28))
datetick('x',6,'keeplimits','keepticks')
text(t1-3,-130,datestr(t1,2))
title(['birthday: ' datestr(birthday,2)])
axis tight
legend('Physical','Emotional','Intellectual')
l2 = legend('Physical','Emotional','Intellectual');
set(l2,'pos',[.07 .60 .18 .10])
