% OSCILLATOR.  Compare four integrators, including classical Runge-Kutta.

F = @(t,y) [y(2); -y(1)];
opts = odeset('reltol',1.e-6,'abstol',1.e-6,'refine',1);
y0 = [1 0];
tspan = [0 2*pi];
h = pi/50;
for k = 1:4
   switch k
      case 1, [t,y] = ode23(F,tspan,y0,opts);
      case 2, [t,y] = ode45(F,tspan,y0,opts);
      case 3, [t,y] = ode113(F,tspan,y0,opts);
      case 4, [t,y] = myrk4(F,tspan,y0,h);
   end
   err(k) = max(abs(y(end,:)-y0));
   cnt(k) = length(t)-1;
end
fprintf('       ode23        ode45       ode113        myrk4\n')
fprintf('%12.2e %12.2e %12.2e %12.2e\n',err) 
fprintf('%12.0f %12.0f %12.0f %12.0f\n',cnt) 
