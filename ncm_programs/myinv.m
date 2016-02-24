function X = myinv(A)
%MYINV  My very own inverse.
%   X = myinv(A) is the inverse of a square matrix A.

[n,n] = size(A);

% Triangular factorization
[L,U,p] = lutx(A);

% Permutation and forward elimination
I = eye(n,n);
X = zeros(n,n);
for k = 1:n
   j = 1:k-1;
   X(k,:) = I(p(k),:) - L(k,j)*X(j,:);
end

% Back substitution
for k = n:-1:1
   j = k+1:n;
   X(k,:) = (X(k,:) - U(k,j)*X(j,:))/U(k,k);
end
