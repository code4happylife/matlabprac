function [L,U,p,q] = lupiv(A)
%LUPIV  Triangular factorization with various pivoting strategies.
%   Ordinarily, [L,U,p] = lupiv(A) produces a unit lower triangular
%   matrix L, an upper triangular matrix U and a permutation vector p,
%   so that L*U = A(p,:)
%
%   With only two output arguments, [L,U] = lupiv(A) does no pivoting.
%   With four output arguments, [L,U,p,q] = lupiv(A) does complete pivoting.

[n,n] = size(A);
p = (1:n)';
q = 1:n;

for k = 1:n-1

   if nargout == 2
      % No pivoting
      i = k;
      j = k;
   elseif nargout <= 3
      % Find index of largest element below diagonal in k-th column
      i = find(abs(A(k:n,k)) == max(abs(A(k:n,k))))+k-1;
      j = k;
   else
      % Find indices of largest element in unreduced matrix.
      [i,j] = find(abs(A(k:n,k:n)) == max(max(abs(A(k:n,k:n)))));
      i = i+k-1;
      j = j+k-1;
   end

   % Skip elimination if pivot is zero
   if (A(i,j) ~= 0)
   
      % Swap pivot row
      if (i ~= k)
         A([k i],:) = A([i k],:);
         p([k i]) = p([i k]);
      end

      % Swap pivot column
      if (j ~= k)
         A(:,[k j]) = A(:,[j k]);
         q([k j]) = q([j k]);
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
L = tril(A,-1) + eye(n,n);
U = triu(A);
