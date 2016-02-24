syms a b c fa fb fc

% Lagrange

x1 = (0-fb)*(0-fc)/((fa-fb)*(fa-fc))*a + ...
        (0-fa)*(0-fc)/((fb-fa)*(fb-fc))*b + ...
        (0-fa)*(0-fb)/((fc-fa)*(fc-fb))*c;
disp('x1 = '), pretty(x1)

% Cramer

x = [a; b; c];
f = [fa; fb; fc];
e = ones(3,1);
x2 =  det([f.^2 f x])/det([f.^2 f e]);
disp('x2 = '), pretty(x2)

% fzerotx, iqi

m = (a - b)/2;
s = fb/fc;
q = fc/fa;
r = fb/fa;
p = s*(2*m*q*(q - r) - (b - c)*(r - 1));
q = (q - 1)*(r - 1)*(s - 1);
p = -p;
d = p/q;
x3 = b+d;
disp('x3 = '), pretty(x3)

simplify(x1-x3)
simplify(x2-x3)
