% TRACKEIGS  Eigenvalue trajectories.
shg
clf
set(gcf,'double','on')
n = 12;
[I,J] = ndgrid(1:n);
axis(10*[-1 1 -1 1])
axis square
hold on
e = 10*ones(n,1);
h = plot(real(e),imag(e),'.','erasemode','none');
t = 0.100;
while any(abs(e) < 11)
   A = 1./(I-J+t);
   e = eig(A);
   set(h,'xdata',real(e),'ydata',imag(e))
   drawnow
   t = t + .005;
end
hold off
