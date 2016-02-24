function x = bslashtx2(A,b)
%BSLASHTX2  Solve linear system (backslash)
%   x = bslashtx2(A,b) solves A*x = b

% Triangular factorization
[L,U,p] = lutx(A);

% Permutation and forward elimination
y = forward(L,b(p));

% Back substitution
x = backsubs(U,y);


% ------------------------------

function x = forward(L,x)
% FORWARD. Forward elimination.
% For lower triangular L, x = forward(L,b) solves L*x = b.
[n,n] = size(L);
for k = 1:n
   x(k) = x(k)/L(k,k);
   i = k+1:n;
   x(i) = x(i) - L(i,k)*x(k);
end

% ------------------------------

function x = backsubs(U,x)
% BACKSUBS.  Back substitution.
% For upper triangular U, x = backsubs(U,b) solves U*x = b.
[n,n] = size(U);
for k = n:-1:1
   x(k) = x(k)/U(k,k);
   i = 1:k-1;
   x(i) = x(i) - U(i,k)*x(k);
end
