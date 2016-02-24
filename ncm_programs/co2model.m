function co2model
% Simulation of the carbon system in atmosphere, shallow, and deep sea.
% See J. C. G. Walker, Numerical adventures with geochemical cycles,
% Oxford University Press, 1991.

% Mass balance equations for:
%   p    = y(1) = atmospheric CO2
%   sigs = y(2) = CO2 in the shallow ocean
%   sigd = y(3) = CO2 in the deep ocean 
%   alks = y(4) = alkalinity in the shallow ocean
%   alkd = y(5) = alkalinity in the deep ocean

% Integration parameters and initial values
tspan = [1000 5000];
y0 = [1.00 2.01 2.23 2.20 2.26];

% Solve ODEs
[t,y] = ode23tx(@co2,tspan,y0);

% CO2 problem, part a, plot output.

format compact
disp('Part a: see figure(1)')
figure(1)
subplot(2,1,2)
plot(t,y(:,1:3))
axis([tspan 0 4.5]);
legend('atmosphere','shallow','deep',4);
xlabel('time (yr)');
ylabel('carbon')
subplot(2,1,1)
plot(t,fuel(t),'k');
ylabel('fuel')
legend('fossil fuel',4)
title('Carbon in the atmosphere and ocean')
drawnow

% CO2 problem, part b, compare initial and final values.
      
disp(' ')
disp('Part b:')
disp('      p        sigs      sigd      alks      alkd')
y0 
yfinal = y(end,:)
ratio = yfinal(1:3)./y0(1:3)
disp('Large increase in the atmosphere.')
disp('Slight increases in the ocean.')

% CO2 problem, part c, maximum atmospheric co2.
% Solve ODEs with gstop and increased accuracy.

opts = odeset('events',@gstop,'reltol',1.e-6);
[tode45,yode45,te,ye] = ode45(@co2,tspan,y0,opts);
disp(' ')
disp('Part c:')
attime = te
maxco2 = ye(1)

% Results:
% attime = 2348
% maxco2 = 4.217

% CO2 problem, part d, observe stiffness.

disp(' ')
disp('Part d: stiffness, see figure(2)')
figure(2)
subplot(2,1,1)
xz = [3000 3100 3100 3000 3000];
yz = [3.12 3.12 3.24 3.24 3.12];
plot(t,y(:,1),'b-',xz,yz,'r-')
text(3150,3.25,'zoom')
subplot(2,1,2)
plot(t,y(:,1))
axis([3000 3100 3.12 3.24])
drawnow

% Problem 7.18, part e, various solvers.

disp(' ')
disp('Part e: compare solvers.')
opts = odeset('reltol',1.e-6,'abstol',1.e-6,'stats','on');
disp(' ')
disp('ode23tx')
tic
[t,y] = ode23tx(@co2,tspan,y0,opts);
ode23txsteps = length(t)
ode23txtime = toc

disp(' ')
disp('ode23')
tic
[t,y] = ode23(@co2,tspan,y0,opts);
ode23time = toc

disp(' ')
disp('ode45')
tic
[t,y] = ode45(@co2,tspan,y0,opts);
ode45time = toc

disp(' ')
disp('ode113')
tic
[t,y] = ode113(@co2,tspan,y0,opts);
ode113time = toc

disp(' ')
disp('ode23s')
tic
[t,y] = ode23s(@co2,tspan,y0,opts);
ode23stime = toc

disp(' ')
disp('ode15s')
tic
[t,y] = ode15s(@co2,tspan,y0,opts);
ode15stime = toc


%------------------------

function f = fuel(t)
% Fossil fuel rates
d = [1000   0.0
     1850   0.0
     1950   1.0
     1980   4.0
     2000   5.0
     2050   8.0
     2080  10.0
     2100  10.5
     2120  10.0
     2150   8.0
     2225   3.5
     2300   2.0
     2500   0.0
     5000   0.0];
f = pchip(d(:,1),d(:,2),t);

%------------------------

function dydt = co2(t,y)
% Mass balance equations for atmospheric CO2, CO2 in the shallow and
% deep ocean, and alkalinity in the shallow and deep ocean.

% Parameters
d = 8.64;
mu1 = 4.95e+2;
mu2 = 4.95e-2;
vs = 0.12;
vd = 1.23;
w = 1.0e-3;
k1 = 2.19e-4;
k2 = 6.12e-5;
k3 = 0.997148;
k4 = 6.79e-2;

% Primary dependent variables
p = y(1);
sigs = y(2);
sigd = y(3);
alks = y(4);
alkd = y(5);

% Equilibrium equations
hco3s = (sigs-sqrt(sigs^2-k3*alks*(2*sigs-alks)))/k3;
co3s = (alks-hco3s)/2;
ps = k4*hco3s^2/co3s;

% Differential equations
f = fuel(t);
dydt = [(ps-p)/d+f/mu1
        ((sigd-sigs)*w-k1-(ps-p)/d*mu2)/vs
        (k1-(sigd-sigs)*w)/vd
        ((alkd-alks)*w-k2)/vs
        (k2-(alkd-alks)*w)/vd];

%------------------------
% gstop function for use in part c.

function [g,isterminal,direction] = gstop(t,y)
% Find zero of derivative of atmospheric co2.
ydot = co2(t,y);
g = ydot(1);
isterminal = 1;
direction = [];

