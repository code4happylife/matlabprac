% TRIPLETS
%    Eigenvalues of the square with multiplicity greater than two.
%    Find integers that can be written as the sum of squares
%    in more than two ways.

for n = 1:325
   m = floor(sqrt(n-1));
   s = n - (1:m).^2;
   k = find(sqrt(s) == round(sqrt(s)));
   if length(k) > 2
      t = sprintf('%3d',n);
      s = s(k);
      for j = 1:length(s)
         t = [t sprintf(' = %d^2+%d^2',sqrt(s(j)),sqrt(n-s(j)))];
      end
      disp(t)
   end
end
