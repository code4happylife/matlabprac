function zdot = longjump(t,z,rho)
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
