% Train whistle

load train
y = y';
Fs 
n = length(y)

t = (1:n)/Fs;
figure(1)
clf reset
plot(t,y)
xlabel('time (seconds)')
ylabel('amplitude')
title('Toot Toot')

Y = fft(y);
Y(1) = [];
f = (1:n/2)/n*Fs;
figure(2)
clf reset
plot(f,abs(Y(1:n/2)));
xlabel('frequency (cycles/second)')
ylabel('power')
title('Train whistle periodogram')

% Pick off six peaks with

% [peaks,q] = ginput(6);

% These six peaks are roughly

format short
peaks = [700 875 1167 2100 2625 3500]'

% The ratios are

format rat
peaks = [700 875 1166.666667 2100 2625 3500]';
ratios = peaks/peaks(1)

disp('peaks(1), peaks(2), peaks(3) are fundamental.')
disp('peaks(4) and peaks(6) are the first two overtones of peaks(1).')
disp('peaks(5) is the first overtone of peaks(2).')
