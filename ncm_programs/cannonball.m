function cannonball

shg
clf
x0 = 0;
y0 = 0;
v0 = 50;
tspan = [0 20];
opts = odeset('events',@events,'refine',6);
n = 17;
T = zeros(n,4);
X = zeros(n,4);
V = zeros(n,4);
LT = zeros(n,4);
w = (0:n-1)'/(n-1);
c = [0*w 2/3*(1-w) 2/3*w];
set(gcf,'defaultaxescolororder',c)

for w = 1:4
   subplot(2,2,w)
   axis([0 200 0 110])
   hold on
   for k = 1:n
      theta0 = k*pi/36;
      z0 = [x0; y0; v0; theta0];
      [t,z] = ode45(@cannonfun,tspan,z0,opts,w);
      plot(z(:,1),z(:,2),'-','color',c(k,:))
      drawnow
      T(k,w) = t(end);
      X(k,w) = z(end,1);
      V(k,w) = z(end,3);
      LT(k,w) = length(t);
   end
   hold off
   switch w
      case 1
         title('None')
      case 2
         title('Head')
      case 3
         title('Tail')
      case 4
         title('Gusty')
   end
end

kays = {};
for k = 1:17
   kays = {kays{:},int2str(5*k)};
end
set(legend(kays),'pos',[.01 .25 .075 .50])

T
X
V 
LT

for w = 1:4
   switch w
      case 1
         fprintf('\nNo wind\n')
      case 2
         fprintf('\nSteady head wind\n')
      case 3
         fprintf('\nIntermittent tail wind\n')
      case 4
         fprintf('\nGusty wind\n')
   end
   k = find(X(:,w) == max(X(:,w)));
   fprintf('  theta0 = %2.0f degrees, t = %6.3f, range = %5.1f, vfinal = %6.3f, nsteps = %2.0f\n', 5*k, T(k,w), X(k,w), V(k,w), LT(k,w));
end

% ---------------------------------------------

function zdot = cannonfun(t,z,wind)
% Cannonball problem.

x = z(1);
y = z(2);
v = z(3);
theta = z(4);

c = 0.2;     % Drag coefficient
rho = 1.29;  % Density of air
s = 0.25;    % Cross section
m = 15;      % Mass
g = 9.8;     % Gravity

switch wind
   case 1
      W = 0;
   case 2
      W = -10;
   case 3
      if rem(floor(t),2) == 0, W = 10; else, W = 0; end
   case 4
      W = 10*randn;
end

xdot = v*cos(theta);
ydot = v*sin(theta);
D = c*rho*s*((xdot-W)^2 + ydot^2)/2;
vdot = -D/m - g*sin(theta);
thetadot = -(g/v)*cos(theta);

zdot = [xdot; ydot; vdot; thetadot];

% ---------------------------------------------

function [value,isterminal,direction] = events(t,y,ignore)
% Locate the time when height passes through zero in
% a decreasing direction and stop integration.  
value = y(2);     % detect height = 0
isterminal = 1;   % stop the integration
direction = -1;   % negative direction
