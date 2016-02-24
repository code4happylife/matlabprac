% Fourier analysis of El Nino data (Southern Oscillation Index)

load elnino.dat
soi = elnino';
n = length(soi)

% Data is monthly from 1962 through 1976.
t = 1962 + (0:n-1)/12;

% Linear trend
c = polyfit(t,soi,1);
trend = polyval(c,t);

% Plot data with trend
figure(1)
plot(t,[soi; trend],'-',t,soi,'k.')
xlabel('year')
ylabel('SOI')
title('Southern Oscillation Index')

% FFT
y = soi - trend;
Y = fft(y);

% Delete first Fourier coefficient, because Y(1) = sum(y) is zero.

Y(1) = [];

% Periodogram with x axis in cycles per month.

figure(2)
pow = abs(Y(1:n/2)).^2;
pmax = 7e4;
f = (1:n/2)/n;
plot([f; f],[0*pow; pow],'c-', f,pow,'b.','linewidth',2,'markersize',16)
axis([0 .5 0 pmax])
xlabel('cycles/month')
ylabel('power')
title('Periodogram')

% Zoom in and invert x axis to years per cycle.

figure(3)
k = 1:16;
pow = pow(k);
mpk = n./k;  % Months per cycle
plot([k; k],[0*pow; pow],'c-',k,pow,'b.','linewidth',2,'markersize',16)
axis([min(k)-1 max(k)+1 0 pmax])
set(gca,'xtick',k)
xticklabels = sprintf('%5.1f|',mpk);
set(gca,'xticklabel',xticklabels)
xlabel('months/cycle')
ylabel('power')
title('Periodogram')

disp('The strongest peak is at 12 months per cycle, i.e. yearly.')
disp('There is another peak spread across three components with 28,')
disp('33.6, and 42 months per cycle, i.e. a little less than three years.')
