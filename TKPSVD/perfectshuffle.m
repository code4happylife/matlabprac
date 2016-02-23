function S = perfectshuffle(p,q)
p=3;
q=3;
% S = perfectshuffle(p,q)
% -------------------------
% Constructs the pq-by-pq perfect shuffle matrix S for a p-by-q matrix A.
%
% S         =   perfect shuffle matrix,
%
% p         =   scalar, number of rows of A,
%
% q         =   scalar, number of columns of A.
% The entries of kron(B,C) and kron(C,B)consists of all possible products
% of a B-matrix entry with a C-matrix entry and this raises the possibility
% that those two Kronecker products are related by a permutation.The
% permutation involved is in fact the perfect shuffle（完全打乱，完美洗牌）.

r=p*q;
I=eye(r);%构造一个r*r的单位矩阵
S=[];
for i=1:q
   S=[S; I(i:q:r,:)];%将I矩阵，按行进行洗牌操作，因为p*q=r
end