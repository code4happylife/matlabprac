% INVEST.  $100 at 6 percent over 10 years.

y0 = 100;
r = .06;
t = 10;

% Yearly

h = 1
y = y0;
for k = 1:t/h
   y = y + h*r*y;
end
yearly = y;

% Monthly

h = 1/12;
y = y0;
for k = 1:t/h
   y = y + h*r*y;
end
monthly = y;

% Midpoint rule

y = y0;
for k = 1:t/h
   s1 = r*y;
   s2 = r*(y + h/2*s1);
   y = y + h*s2;
end
midpoint = y;

% Trapezoid rule

y = y0;
for k = 1:t/h
   s1 = r*y;
   s2 = r*(y + h*s1);
   y = y + h*(s1+s2)/2;
end
trapezoid = y;

% BS23

y = y0;
s4 = r*y;
for k = 1:t/h
   s1 = s4;
   s2 = r*(y + h*s1/2);
   s3 = r*(y + 3*h*s2/4);
   y = y + h/9*(2*s1+3*s2+4*s3);
   s4 = r*y;
end
bs23 = y;

% Compound continuously

continuous = exp(r*t)*y0;

[yearly; monthly; midpoint; trapezoid; bs23; continuous]
