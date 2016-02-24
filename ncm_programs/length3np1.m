function L = length3np1(n)
% LENGTH3NP1  Length of 3n+1 sequence.
% LENGTH3NP1(N)
% Does not save the sequence, so uses little memory.
L = 1;
while n > 1
   if rem(n,2)==0
      n = n/2;
   else
      n = 3*n+1;
   end
   L = L + 1;
end
