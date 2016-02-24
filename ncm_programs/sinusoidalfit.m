% Problem 5_6.  Least squares fit with outliner and sinusoid.

t = 1:25;
y = [ 5.0291  6.5099  5.3666  4.1272  4.2948 ...
      6.1261 12.5140 10.0502  9.1614  7.5677 ...
      7.2920 10.0357 11.0708 13.4045 12.8415 ...
     11.9666 11.0765 11.7774 14.5701 17.0440 ...
     17.0398 15.9069 15.4850 15.5112 17.6572];
figure(1)
plot(t,y,'o');

t = t';
y = y';
A = [ones(size(t)) t];
c = A\y
f = A*c;
figure(2)
bar(t,y-f);
title('(t(7),y(7)) is the outlier')

t7 = t(7);
y7 = y(7);
t(7) = [];
y(7) = [];
A = [ones(size(t)) t sin(t)];
c = A\y
s = (0:.1:26)';
f = [ones(size(s)) s sin(s)]*c;
figure(3)
plot(t,y,'o',s,f,'-',t7,y7,'*')
axis tight
