% Various least squares fits to erf(t)

m = 10;
nmax = 10;
tmax = 1;
t = (0:tmax/m:tmax)';
mt = length(t);
y = erf(t);
u = (0:tmax/(20*m):tmax)';
mu = length(u);

% Polynomial fit

err = zeros(4,nmax);
for n = 1:nmax
   c = polyfit(t,y,n);
   v = polyval(c,u);
   err(1,n) = max(abs(v-erf(u)));
end

% Odd degree polynomial fit

for n = 1:nmax
   A = zeros(mt,n);
   A(:,1) = t;
   for k = 2:n
      A(:,k) = t.^2 .* A(:,k-1);
   end
   c = A\y;
   w = u;
   z = u.^2;
   v = c(1)*u;
   for k = 2:n
      w = z.*w;
      v = v + c(k).*w;
   end
   err(2,n) = max(abs(v-erf(u)));
end

% erf(x) ~= c(1) + exp(-x^2)*(c(2) + c(3)/(1+t) + ... + c(n)/(1+t)^(n-2))

for n = 2:nmax
   A = zeros(mt,n);
   A(:,1) = 1;
   A(:,2) = exp(-t.^2);
   z = 1./(1+t);
   for k = 3:n
      A(:,k) = A(:,k-1).*z;
   end
   c = A\y;
   w = exp(-u.^2);
   v = c(1) + c(2)*w;
   z = 1./(1+u);
   for k = 3:n
      w = z.*w;
      v = v + c(k).*w;
   end
   err(3,n) = max(abs(v-erf(u)));
end

% Polynomial fit to expcx(t) = exp(t^2)*(1-erf(t))

z = erfcx(t);
for n = 1:nmax
   c = polyfit(t,z,n);
   v = 1-exp(-u.^2).*polyval(c,u);
   err(4,n) = max(abs(v-erf(u)));
end
   
clf
shg
n = 1:nmax;
semilogy(n,err,'o')
axis([0 nmax+1 1.e-12 1])
legend('poly','odd','exp(-t^2)','erfcx(t)',3)
title(['tmax = ' num2str(tmax)])
