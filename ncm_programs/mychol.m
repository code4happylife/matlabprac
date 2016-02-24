function R = mychol(A)
%MYCHOL My own version of Cholesky factorization.
%   R = mychol(A), for a positive definite matrix A,
%   is an upper triangular matrix R so that R'*R = A.

if ~isequal(A',A)
   error('Matrix must be symmetric or Hermitian')
end
[n,n] = size(A);
R = zeros(n,n);
for k = 1:n
   i = (1:k-1)';
   j = k+1:n;
   s = A(k,k) - R(i,k)'*R(i,k);
   if s > 0
      R(k,k) = sqrt(s);
   else
      error('Matrix must be positive definite')
   end
   R(k,j) = (A(k,j) - R(i,k)'*R(i,j))/R(k,k);
end
