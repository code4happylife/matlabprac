% BEAMON  Bob Beamon's long jump

function longjump

format compact
disp('    v0       theta0     rho       distance')
clf

% Initial angle 22.5 degrees throughout.

theta0 = pi/8;

% Nominal jump at high altitude

rho = .94;
v0 = 10;
xf = f1(v0,theta0,rho);
disp([v0 theta0/pi*180 rho xf])
p1(v0,theta0,rho,[0 2/3 0])

% Air density at sea level

rho = 1.29;
xf = f1(v0,theta0,rho);
disp([v0 theta0/pi*180 rho xf])
p1(v0,theta0,rho,[0 0 0])

% Choose initial velocity so that distance is Beamon's record.

rho = .94;
xf = 8.90;
v0 = fzerotx(@f2,[9,12],theta0,rho,xf);
xf = f1(v0,theta0,rho);
disp([v0 theta0/pi*180 rho xf])
p1(v0,theta0,rho,[0 3/4 3/4])

% Beamon at sea level

rho = 1.29;
xf = f1(v0,theta0,rho);
disp([v0 theta0/pi*180 rho xf])
p1(v0,theta0,rho,[0 0 3/4])

% Zoom box
subplot(2,1,1)
plot([6.2 9.2 9.2 6.2 6.2],[-.04 -.04 .39 .39 -.04],'k:')
hold off

% ------ ode ------

function zdot = longjumpode(t,z,rho)
% Long jump ode

c = 0.72;    % Drag coefficient
s = 0.5;     % Cross section
m = 80;      % Mass
g = 9.81;    % Gravity

x = z(1);
y = z(2);
v = z(3);
theta = z(4);
xdot = v*cos(theta);
ydot = v*sin(theta);
D = c*rho*s*(xdot^2 + ydot^2)/2;
vdot = -D/m - g*sin(theta);
thetadot = -(g/v)*cos(theta);
zdot = [xdot; ydot; vdot; thetadot];

% ------ f1.  Distance function.  How long is the jump? ------

function xf = f1(v0,theta0,rho);
x0 = 0;
y0 = 0;
z0 = [x0; y0; v0; theta0];
opts = odeset('events',@g1);
tspan = [0 2];
[t,z,tf,zf] = ode45(@longjumpode,tspan,z0,opts,rho);
xf = zf(1);

% ------ g1.  Gstop function.  Stop when height = 0. ------

function [gstop,isterm,dir] = g1(t,y,ignore)
gstop = y(2);
isterm = 1;
dir = -1;

% ------ f2.  Compare distance to specified distance. ------

function miss = f2(v0,theta0,rho,xf);
miss = f1(v0,theta0,rho) - xf;

% ------ p1.  Plot jump ------

function p1(v0,theta0,rho,color);
x0 = 0;
y0 = 0;
z0 = [x0; y0; v0; theta0];
opts = odeset('events',@g1);
tspan = 0:.02:2;
[t,z] = ode45(@longjumpode,tspan,z0,opts,rho);

subplot(2,1,1)
h = plot(z(:,1),z(:,2));
axis([0 10 -.25 1.5])
set(h,'color',color)
hold on

subplot(2,1,2)
h = plot(z(:,1),z(:,2));
axis([6.2 9.2 -.04 .39])
set(h,'color',color)
hold on
