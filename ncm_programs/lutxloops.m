function [L,U,p] = lutxloops(A)
%LU Triangular factorization
%   [L,U,p] = lup(A) produces a unit lower triangular matrix L,
%   an upper triangular matrix U and a permutation vector p,
%   so that L*U = A(p,:)

[n,n] = size(A);
p = (1:n)';

for k = 1:n-1

   % Find index of largest element below diagonal in k-th column
   m = k;
   for i = k+1:n
      if abs(A(i,k)) > abs(A(m,k))
         m = i;
      end
   end

   % Skip elimination if column is zero
   if (A(m,k) ~= 0)
   
      % Swap pivot row
      if (m ~= k)
         for j = 1:n;
            A([k m],j) = A([m k],j);
         end
         p([k m]) = p([m k]);
      end

      % Compute multipliers
      for i = k+1:n;
         A(i,k) = A(i,k)/A(k,k);
      end

      % Update the remainder of the matrix
      for j = k+1:n;
         for i = k+1:n;
            A(i,j) = A(i,j) - A(i,k)*A(k,j); 
         end
      end
   end
end
    
% Separate result
L = tril(A,-1) + eye(n,n);
U = triu(A);
