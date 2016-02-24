% BELLSHAPE  Solve d2u/dx2 = exp(-x^2) on [-1,1]

% Numerically

shg
clf
n = 63;
h = 2/(n+1);
x = (-1:h:1)';
f = exp(-x.^2);
e = ones(n,1);
u = zeros(n+2,1);
A = spdiags([e -2*e e],[-1 0 1],n,n)/h^2;
j = 2:n+1;
u(j) = A\f(j);
plot(x,u);
residual = max(abs(diff(u,2)/h^2 - exp(-x(j).^2)))

% Symbolically

U = dsolve('D2u = exp(-t^2)','u(-1) = 0','u(1) = 0');
U = simple(U);
pretty(U)

% Compare

error = max(abs(subs(U,x) - u))

