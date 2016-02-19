A=[1 3 5;7 9 11;2 4 6]
B=[1 2 4;5 8 10;1 2 5]
%definition of matrix-vector multiplication
C=zeros(3,3)
for i=1:3
    for j=1:3
        for s=1:3
            C(i,j)=C(i,j)+A(i,s)*B(s,j)
        end
    end
end
%inner product version of matrix multiplication
C=zeros(3,3)
for i=1:3
    for j=1:3
        C(i,j)=A(i,1:3)*B(1:3,j)
    end
end
% It is immediately seen that the the loop variables can be permuted in 3! = 6 di?erent
% ways, and we can write a generic matrix multiplication code:
% for ...
% for ...
% for ...
% C(i,j)=C(i,j)+A(i,s)*B(s,j)
% end
% end
% end
% The matrix A is accessed by columns and B by scalars. 
% This access pattern can be illustrated as
C=zeros(3,3)
for j=1:3
    for s=1:3
        C(1:3,j)=C(1:3,j)+A(1:3,s)*B(s,j)
    end
end
% In another permutation we let the s-loop be the outermost:
% This can be illustrated as follows. Let a¡¤ k denote the column vectors of A and let
% b T
% k ¡¤
% denote the row vectors of B. Then matrix multiplication can be written as
% follow
% This is the outer product form of matrix multiplication. Remember that the outer
% product follows the standard de?nition of matrix
% multiplication:C(1:3,j) is matrix rather than column vectors.
C=zeros(3,3)
for s=1:3
    for j=1:3
        C(1:3,j)=C(1:3,j)+A(1:3,s)*B(s,j)
    end
end



