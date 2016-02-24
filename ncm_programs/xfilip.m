% Filip data set

load filip.dat
y = filip(:,1);
x = filip(:,2);
n = length(x);
d = 10;
p = d+1;
u = (-8.78:.01:-3.13)';

% NIST reference coefficients

beta0 = flipud([ ...
        -1467.48961422980
        -2772.17959193342
        -2316.37108160893
        -1127.97394098372
        -354.478233703349
        -75.1242017393757
        -10.8753180355343
        -1.06221498588947
        -0.670191154593408E-01
        -0.246781078275479E-02
        -0.402962525080404E-04]);
p0 = polyval(beta0,u);
f0 = polyval(beta0,x);
r0 = norm(y-f0);

% Polyfit

beta1 = polyfit(x,y,d)';
f1 = polyval(beta1,x);
r1 = norm(y-f1);
rsd0 = 0.334801051324544E-02
rsd1 = r1/sqrt(n-p)
p1 = polyval(beta1,u);

% Vandermonde and backslash

X = vander(x);
X = X(:,n-d:n);
beta2 = X\y;
f2 = polyval(beta2,x);
r2 = norm(y-f2);
p2 = polyval(beta2,u);

% Pseudoinverse

beta3 = pinv(X)*y;
f3 = polyval(beta3,x);
r3 = norm(y-f3);
p3 = polyval(beta3,u);

% Normal equations (reject)

beta4 = inv(X'*X)*X'*y;
f4 = polyval(beta4,x);
r4 = norm(y-f4);
p4 = polyval(beta4,u);

% Centered

mu = mean(x);
sigma = std(x);
t = (x-mu)/sigma;
gamma5 = polyfit(t,y,d)';
f5 = polyval(gamma5,t);
r5 = norm(y-f5);
p5 = polyval(gamma5,(u-mu)/sigma);
beta5 = sym2poly(expand(((sym('x')-mu)/sigma).^(d:-1:0)*gamma5)).';

format long e
disp('beta = ')
disp('    NIST                      Polyfit                   Centered')
disp([beta0 beta1 beta5])
normr = [r0 r1 r5]

disp('beta = ')
disp('    NIST                      Vandermonde               Pseudoinverse')
disp([beta0 beta2 beta3])
normr = [r0 r2 r3]

clf
shg
plot(x,y,'.',u,[p0 p2])
set(legend('Data','Rank 11','Rank 10',4),'pos',[.6 .3 .18 .11])
title('NIST Filip data set')
% print -depsc2 filip.eps
