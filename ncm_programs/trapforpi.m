% trapforpi.m  Solution to problem 6.3
a = -1;
b = 1;
f = @(x) 2./(1+x.^2);
fprintf('  n      h             T(h)')
fprintf('           e = T(h)-pi       e/h^2\n')
for n = [11 21 41 81 161 321]
   h = (b-a)/(n-1);
   x = a + h*(1:n-2);
   T = h/2*f(a) + h*sum(f(x)) + h/2*f(b);
   e = T - pi;
   fprintf('%3d %9.5f %19.14f %15.4e %15.9f\n',n,h,T,e,e/h^2')
end
