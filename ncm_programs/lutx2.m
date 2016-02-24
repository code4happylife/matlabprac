function [L,U,p] = lutx(A)
%LUTX2  Triangular factorization, textbook version.
%   With three output arguments, [L,U,p] = LUTX2(A) produces
%   a unit lower triangular matrix L, an upper triangular
%   matrix U, and a permutation vector p, so that L*U = A(p,:)
%
%   With two output arguments, [L,U] = LUTX2(A) produces a
%   "psychologically lower triangular matrix" L (i.e. a product
%   of lower triangular and permutation matrices), and an
%   upper triangular U so that L*U = A.
%
%   With one output argument, LUTX2(A) returns a single matrix
%   containing L-I+U.  The pivot information is lost.

[n,n] = size(A);
p = (1:n)';

for k = 1:n-1

   % Find index of largest element below diagonal in k-th column
   [r,m] = max(abs(A(k:n,k)));
   m = m+k-1;

   % Skip elimination if column is zero
   if (A(m,k) ~= 0)
   
      % Swap pivot row
      if (m ~= k)
         A([k m],:) = A([m k],:);
         p([k m]) = p([m k]);
      end

      % Compute multipliers
      i = k+1:n;
      A(i,k) = A(i,k)/A(k,k);

      % Update the remainder of the matrix
      j = k+1:n;
      A(i,j) = A(i,j) - A(i,k)*A(k,j); 
   end
end

% Separate result
if nargout == 3
   L = tril(A,-1) + eye(n,n);
elseif nargout == 2
   L(p,:) = tril(A,-1) + eye(n,n);
else
   L = A;
end
U = triu(A);
