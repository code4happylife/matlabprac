function ferntrajs
% FERNTRAJS.  Fern trajectories.

x0 = [-1 5]';
n = 30000;
markersize = 1;

x = [.5; .5];
xs = zeros(2,n);
xs(:,1) = x;
p  = [ .85  .92  .99  1.00];
A1 = [ .85  .04; -.04  .85];  b1 = [0; 1.6];
A2 = [ .20 -.26;  .23  .22];  b2 = [0; 1.6];
A3 = [-.15  .28;  .26  .24];  b3 = [0; .44];
A4 = [  0    0 ;   0   .16];
for j = 2:n
   r = rand;
   if r < p(1)
      x = A1*x + b1;
   elseif r < p(2)
      x = A2*x + b2;
   elseif r < p(3)
      x = A3*x + b3;
   else
      x = A4*x;
   end
   xs(:,j) = x;
end

shg
set(gcf,'color','white')
darkgreen = [0 .5 0];
plot(xs(1,:),xs(2,:),'.','markersize',markersize,'color',darkgreen);
axis([-3 3 0 10])
axis off

x0 = x;
hold on
I = eye(2,2);
for k = 1:4
   if k == 1, n = 20; else, n = 5; end
   x = x0;
   xs = zeros(2,n);
   xs(:,1) = x;
   for j = 2:n
      switch k
         case 1
            x = A1*x + b1;
         case 2
            x = A2*x + b2;
         case 3
            x = A3*x + b3;
         case 4
            x = A4*x;
      end
      xs(:,j) = x;
   end
   switch k
      case 1
         z = (I - A1)\b1;
         dz = [.15; 0];
      case 2
         z = (I - A2)\b2;
         dz = [0; -.35];
      case 3
         z = (I - A3)\b3;
         dz = [-.35; 0];
      case 4
         z = [0; 0];
         dz = [-.25; 0];
   end
   xs(:,n+1) = z;
   plot(xs(1,:),xs(2,:),'b-o')
   z = z + dz;
   text(z(1),z(2),num2str(k))
end
plot(x0(1),x0(2),'r*')
hold off
