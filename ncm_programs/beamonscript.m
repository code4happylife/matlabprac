rho = 1.29;
x0 = 0;
y0 = 0;
theta0 = pi/8;
z0 = [x0; y0; v0; theta0];
opts = odeset('events',@gstop);
tspan = [0 2];
[t,z,tf,zf] = ode45(@longjump,tspan,z0,opts,rho);
xf = zf(1)

