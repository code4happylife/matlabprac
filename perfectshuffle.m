function S = perfectshuffle(p,q)
r=p*q;
I=eye(r);
S=[];
for i=1:q
   S=[S; I(i:q:r,:)];%i经过特定步长的增加，所得到的可以小于或等于r
end


% S = perfectshuffle(p,q)
% -------------------------
% Constructs the pq-by-pq perfect shuffle matrix S for a p-by-q matrix A.
%
% S         =   perfect shuffle matrix,
%
% p         =   scalar, number of rows of A,
%
% q         =   scalar, number of columns of A.
% There exist a matrix and do some shuffle to the matrix.
% Mess it up and therefore get a "random" version of the original 
% diagonal matrix.
% It uses for loop to get S matrix.
% We should pay more attention on the [] operation.




