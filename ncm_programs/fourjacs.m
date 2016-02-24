% Fourjacs
% Four different equations with solution sin(x).
% Jacobian for F2 is singular at pi/2.

for tf = [pi/2 pi]

   figure

   opts = odeset('abstol',1.e-9,'reltol',1.e-9,'stats','on');
   tspan = [0 tf];
   
   subplot(2,2,1)
   F1 = @(t,y) cos(t)
   ode45(F1,tspan,0,opts);
   axis([tspan 0 1.1])
   
   subplot(2,2,2)
%  if tf == pi/2
      F2 = @(t,y) sqrt(1-y.^2)
%  else
%     F2 = @(t,y) sqrt(abs(1-y.^2))
%  end
   ode45(F2,tspan,0,opts);
   axis([tspan 0 1.1])
   
   subplot(2,2,3)
   F3 = @(t,y) [y(2); -y(1)]
   ode45(F3,tspan,[0,1],opts);
   axis([tspan 0 1.1])
   
   subplot(2,2,4)
   F4 = @(t,y) [y(2); -sin(t)]
   ode45(F4,tspan,[0,1],opts);
   axis([tspan 0 1.1])
end

