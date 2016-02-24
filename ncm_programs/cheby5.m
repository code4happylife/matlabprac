%CHEBY5  Different ways to compute polynomials.
%  Here are several different ways to compute the fifth Chebyshev
%  polynomial, both symbolically and numerically.

% Symbolic manipulation

syms x
phi = (1+sqrt(sym(5)))/2

% Monomial basis
T = 16*x^5-20*x^3+5*x;

% Horner
T(2,1) = ((((16*x + 0)*x - 20)*x + 0)*x + 5)*x + 0;

% Lagrangian interpolation
s = [-1 -phi/2 -(phi-1)/2 (phi-1)/2 phi/2 1];
y = [-1 1 -1 1 -1 1];
T(3) = polyinterp(s,y,x);

% Roots
z = [-sqrt((2+phi)/4) -sqrt((3-phi)/4) 0 ...
      sqrt((3-phi)/4) sqrt((2+phi)/4)];
T(4) = 16*prod(x-z);

% Chebyshev definition
T(5) = expand(cos(5*acos(x)));

% Chebyshev recurrence
C = [1; x];
for n = 2:5
   C(n+1) = 2*x*C(n) - C(n-1);
end
T(6) = C(6);

% Maple knows about orthogonal polynomials
T(7) = maple('orthopoly[T]',5,x);

% Simplify to see if all representations are the same.
T = simple(simple(T))

% ------------------------------------

% Numerical evaluation

x = -1:1/128:1;
phi = (1+sqrt(5))/2

% Monomial basis
T = 16*x.^5-20*x.^3+5*x;

% Horner
T(2,:) = ((((16*x + 0).*x - 20).*x + 0).*x + 5).*x + 0;

% Lagrangian interpolation
s = [-1 -phi/2 -(phi-1)/2 (phi-1)/2 phi/2 1];
y = [-1 1 -1 1 -1 1];
T(3,:) = polyinterp(s,y,x);

% Roots
z = [-sqrt((2+phi)/4) -sqrt((3-phi)/4) 0 ...
      sqrt((3-phi)/4) sqrt((2+phi)/4)];
[Z,X] = ndgrid(z,x);
T(4,:) = 16*prod(X-Z);

% Chebyshev definition
T(5,:) = cos(5*acos(x));

% Chebyshev recurrence
C = [ones(size(x)); x];
for n = 2:5
   C(n+1,:) = 2*x.*C(n,:) - C(n-1,:);
end
T(6,:) = C(6,:);

% Check if all the representations are equal to within roundoff error
err = max(max(abs(diff(T))))

% Plot results
plot(x,T(1,:),'-',s,y,'o',z,0*z,'o')
axis(1.1*axis)
