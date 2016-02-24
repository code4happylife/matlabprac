% Trefethen 100-digit challenge, problem #1.

f = @(x) cos(log(x)./x)./x;
g = @(x,k) log(x)./x+(k-1/2)*pi';

% Compute T = int(f(x),0,1), but can't do that in one step.
% Use g(x) to find x_k = solution to log(x)/x = -(k-1/2)*pi.
% Then f(x_k) = 0.
% Compute Q_k = int(f(x),x_k,x_{k-1}).
% Then T = sum(Q_k).
% The Q_k's alternate in sign, so the partial sums alternate
% between above and below T.  The infinite sum converges slowly
% so use Aitken's del-squared acceleration.
   
format long
format compact
tol = 1.e-15;
opts = optimset;
havelambertw = exist('lambertw')==2;

k = 1;
xkm1 = 1;
if havelambertw
   p = pi/2;
   xk = lambertw(p)/p;
else
   xk = fzero(g,[.4 .5],opts,1);
end
Qk = quadl(f,xk,xkm1,tol);

splus = Qk;  % sminus < T < splus
sminus = 0;
tkm2 = 0;
tkm1 = 0;
tk = Qk;

for k = 2:10000
   xkm1 = xk;
   if havelambertw
      p = (2*k-1)*pi/2;
      xk = lambertw(p)/p;
   else
      xk = fzero(g,[xk/2 xk],opts,k);
   end

   Qk = quadl(f,xk,xkm1,tol);
   if Qk > 0 
      splus = sminus + Qk;
   else
      sminus = splus + Qk;
   end
   t = (sminus+splus)/2;
   tkm2 = tkm1;
   tkm1 = tk;
   tk = t;
   tkbar = tkm2-(tk-tkm1)^2/(tk-2*tkm1+tkm2);
   if mod(k,100) == 0
      fprintf(1,'%6.0d %15.12f %15.12f %15.12f %15.12f\n', ...
                 k,sminus,tk,tkbar,splus) 
   end
end
