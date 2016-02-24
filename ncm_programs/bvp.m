% Boundary value problem

function Y = bvpmain

% (a)  Shooting

disp('Shooting')
subplot(3,2,1)
format long
eta = fzerotx(@bvpfun1,[1 2])
[x,y] = ode45(@bvpode,[0:.01:1],[0; eta]);
plot(x,y(:,1))
title('Shooting Method')
Y = y(:,1);

% (b) Quadrature

disp('Quadrature')
subplot(3,2,2)
kappa = fzerotx(@bvpfun2g,[.75 1])
eta = sqrt(2*kappa);
[x,y] = ode45(@bvpode,[0:.01:1],[0; eta]);
plot(x,y(:,1))
title('Quadrature')
Y(:,2) = y(:,1);

% (c)  Finite differences, linear iteration

disp('Linear iteration')
subplot(3,2,3)
n = 99
h = 1/(n+1);
e = ones(n,1);
A = full(spdiags([e -2*e e],[-1 0 1],n,n));
b = zeros(n,1);
b(n) = 1;
y =(1:n)'*h;
z = zeros(n,1);
iter = 0;
while h*norm(y-z) > 1.e-14
   z = y;
   y = A\(h^2*(y.^2-1) - b);
   iter = iter+1;
end
iter
x =(0:n+1)'*h;
y = [0; y; 1];
plot(x,y)
title('Linear iteration')
Y(:,3) = y;


% (d)  Finite differences, Newton's method

disp('Newton''s method')
subplot(3,2,4)
n = 99
h = 1/(n+1);
e = ones(n,1);
A = spdiags([e -2*e e],[-1 0 1],n,n);
b = zeros(n,1);
b(n) = 1;
y = (1:n)'*h;
z = ones(n,1);
iter = 0;
while h*norm(y-z) > 1.e-14
   z = y;
   f = A*y + b - h^2*(y.^2 - 1);
   J = A - h^2*spdiags(2*y,0,n,n);
   y = y - J\f;
   iter = iter+1;
end
iter
x =(0:n+1)'*h;
y = [0; y; 1];
plot(x,y)
title('Newton''s Method')
Y(:,4) = y;


% (e)  MATLAB boundary value problem solver

disp('bvp4c')
subplot(3,2,5)
n = 99;
h = 1/(n+1);
x = (0:n+1)*h;
s.x = x;
s.y = [x; ones(size(x))];
s = bvp4c(@bvpode,@bvpbc,s);
x = s.x;
y = s.y;
plot(x,y(1,:))
title('BVP4C')
Y(:,5) = y(1,:)';


%-------------------------------------------

function ydot = bvpode(x,y)
ydot = [y(2); y(1)^2-1];

%-------------------------------------------

function y1 = bvpfun1(eta)
% Boundary value problem, shooting method.
% bvpfun1(eta) solves
%   y'' = y^2 - 1
%   y(0) = 0
%   y'(0) = eta
% and returns y(1)-1

[x,y] = ode45(@bvpode,[0 1],[0; eta]);
y1 = y(end,1)-1;

%-------------------------------------------

function g = bvpfun2g(kappa)
% Boundary value problem, quadrature method.
% g = bvpfun2(kappa) evaluates the integral
%   g(kappa) = \int_0^1 1/sqrt(2*(kappa+y^3/3-y)) dy
% and returns g-1.

g = quadtx(@bvpfun2h,0,1,1.e-8,kappa) - 1;

%-------------------------------------------

function h = bvpfun2h(y,kappa)
% Boundary value problem, quadrature method.
% h = bvpfun3(kappa) evaluates the integrand
%   h(kappa) = 1/sqrt(2*(kappa+y^3/3-y));

h = 1/sqrt(2*(kappa+y^3/3-y));

%-------------------------------------------

function bc = bvpbc(y0,y1)
bc = [y0(1); y1(1)-1];
