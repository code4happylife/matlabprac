function [S,Q,R,I,Ut,j]=ttr12tucker(U,sigmas,V,n)
% [S,Q]=ttr12tucker(U,sigmas,V,n)
% -------------------------------
% Converts a tensor A of size n in the TTr1SVD format to the Tucker (HOSVD) format.
% Usually results in a more sparse core S compared to traditional methods
% (e.g. Alternating Least Squares). By setting entries of the sigmas vector
% to zero the Tucker format of a truncated TTr1 will be obtained.
%
% S         =   d-way array, the core tensor in the Tucker decomposition,
%
% Q         =   cell, each Q{k} contains an orthogonal matrix of the
%               kth mode outer product vectors, 
%
% U         =   cell, contains the U vectors of each of the SVDs in the
%               TTr1 tree,张量的n模展开矩阵进行SVD分解得到的左矩阵
%
% sigmas    =   vector, contains the final singular values in the linear
%               combination of rank-1 terms.
%
% V         =   cell, contains the V vectors of each of the SVDs in the
%               TTr1 tree,
%
% n         =   vector, size of the original tensor A.


d=length(n);
indices=leave2ind(find(sigmas),n); % only handle nonzero sigmas
indices
nonzerosigmas=sigmas(find(sigmas));

Ut=cell(1,d);
% Concatenate all U and V vectors along each mode
for i=d:-1:2
    for j=min(indices(:,2*(i-2)+1)):max(indices(:,2*(i-2)+1))%注意遍历的是同一列的
        I=find(indices(:,2*(i-2)+1)==j);
        %增加一列数据
        Ut{d-i+1}=[Ut{d-i+1} U{j}(:,intersect(indices(I,2*(i-2)+2),indices(I,2*(i-2)+2)))];
        if i==2
            Ut{d}=[Ut{d} V{j}(:,intersect(indices(I,2*(i-2)+2),indices(I,2*(i-2)+2)))];
        end
    end
end

for i=1:d,
    [Q{i},R{i}]=qr(Ut{i},0); % economic QR factorizations of the concatenated U's and V's
    r(i)=rank(Ut{i});
    % distribute the R vectors in a new R cell according to the indices
    %% in order for the summation over the nonzero sigmas to be correct
    temp=R{i};
    R{i}=zeros(size(temp,1),length(nonzerosigmas));
    if i==d
        counter=1;
        Rveccounter=1;
        while counter<=size(indices,1)
            I=intersect(find(indices(:,1)==indices(counter,1)),find(indices(:,2)==indices(counter,2)));
            counter=I(end)+1; 
            R{i}(:,I)=temp(:,Rveccounter)*ones(1,length(I));            
            Rveccounter=Rveccounter+1;
        end
    else
        counter=1;
        Rveccounter=1;
        while counter<=size(indices,1)
            I=intersect(find(indices(:,end-(i-1)*2-1)==indices(counter,end-(i-1)*2-1)),find(indices(:,end-(i-1)*2)==indices(counter,end-(i-1)*2)));
            counter=I(end)+1; 
            R{i}(:,I)=temp(:,Rveccounter)*ones(1,length(I));            
            Rveccounter=Rveccounter+1;
        end
    end
end

% construct the core tensor
S=zeros(prod(r),1);
for k=1:length(nonzerosigmas)
    temp=nonzerosigmas(k)*kron(R{d}(:,k),R{d-1}(:,k));
    for i=d-2:-1:1
        temp=kron(temp,R{i}(:,k));
    end
    S=S+temp;    
end
tol=max(n)*eps(max(nonzerosigmas));
S(abs(S)<tol)=0;
S=reshape(S,r);
end


% i=3;                 i=2; j=(2:4)%j依次等于2,3,4进行遍历I={1，2};{3，4};{5，6}
% 2*(i-2)+1 = 3;        2*(i-2)+1 = 1;      %内层循环针对j，j=3,I={3,4}
% 2*(i-2)+2 = 4;        2*(i-2)+2 = 2;
% j=1;
% d-i+1=1;			d-i+1=2;
% %Ut{1}里面添加U{1}的各个前三列的元素
% for i=d:-1:2
%     for j=min(indices(:,2*(i-2)+1)):max(indices(:,2*(i-2)+1))%注意遍历的是同一列的；j=min(indices(:,2*(i-2)+1)):max(indices(:,2*(i-2)+1))对倒数第二列进行比较，指明了是针对第几次svd运算
%         I=find(indices(:,2*(i-2)+1)==j);%找到针对这些svd的行索引，相当于叶子节点的索引
%         %增加一列数据;I 返回的是索引;intersect返回的是元素                                 %Ut{1}=U{1}(:,1,2,3)指代第一次分解得到的U{1}矩阵，后面才开始存储V矩阵是因为，后面的svd分解，其实都是由第一次分解得到的V{1}矩阵进行svd分解得到的结果
%         Ut{d-i+1}=[Ut{d-i+1} U{j}(:,intersect(indices(I,2*(i-2)+2),indices(I,2*(i-2)+2)))];% U{j}(:,intersect(indices(I,2*(i-2)+2),indices(I,2*(i-2)+2)))，
%         if i==2																			   %intersect求交集，对于indices矩阵的最后一列，其交集，降低重复，indices矩阵本身就是索引，那么对于求得的交集，仍然是索引
%             Ut{d}=[Ut{d} V{j}(:,intersect(indices(I,2*(i-2)+2),indices(I,2*(i-2)+2)))];	   %Ut{2}=U{2}(:1，2)指代的是整个的U{2}矩阵，索引的列对应的非常完美	Ut{2}=[Ut{2} U{3}(:1，2)]	Ut{2}=[Ut{2} U{4}(:1，2)]
%         end                                                                                %Ut{3}=V{2}(:1，2)指代的是整个的V{2}矩阵,Ut{3}=[Ut{3} V{3}(:1,2)],Ut{3}=[Ut{3} V{4}(:1,2)]
%     end
% end
% /*-0.623182942157870	0.777811415583596	0.0815623834331965
% -0.539833734921171	-0.352350738886366	-0.764479231534311
% -0.565882107195943	-0.520440542803954	0.639467811669521*/
% Ut{2}=U{2}
% Ut{2}=[ Ut{2}U{3} ]
% Ut{2}=[ Ut{2}U{4} ]
% Ut{3}=V{2}(:)
% Ut{3}=[ Ut{3}V{3} ]
% Ut{3}=[ Ut{3}V{4} ]
% 
% find (indices(:,end-1)==indices(1,end-1))
% find (indices(:,end)==indices(1,end))
% while counter<=size(indices,1)
%             I=intersect(find(indices(:,end-(i-1)*2-1)==indices(counter,end-(i-1)*2-1)),find(indices(:,end-(i-1)*2)==indices(counter,end-(i-1)*2)));
%             counter=I(end)+1; 
%             R{i}(:,I)=temp(:,Rveccounter)*ones(1,length(I));            
%             Rveccounter=Rveccounter+1;
%         end
% 
% find (indices(:,end-3)==indices(1,end-3))
% find (indices(:,end-2)==indices(1,end-2))
% 
% 
% for i=1:d,
%     [Q{i},R{i}]=qr(Ut{i},0); % economic QR factorizations of the concatenated U's and V's Ut{1},Ut{2},Ut{3}进行qr分解
%     r(i)=rank(Ut{i});%Ut的rank依次是[3,4,2]
%     % distribute the R vectors in a new R cell according to the indices
%     %% in order for the summation over the nonzero sigmas to be correct
%     temp=R{i};%每次循环将R{1}R{2}R{3}备份到temp临时变量中
%     R{i}=zeros(size(temp,1),length(nonzerosigmas));%对新的R{1}R{2}R{3}进行初始化操作；R{1}R{2}R{3}依次是(3*6)(4*6)(2*6)的矩阵
%     if i==d                                        %i=3,最后执行，将整个Ut{3}QR分解之后的结果的R{3}放到R{3}里面
%     	%%-1	1.11022302462516e-16	-0.243845979509937	0.969813970963937	0.0238799428886511	-0.999714833503852
%         %%  0	1	-0.969813970963937	-0.243845979509937	0.999714833503852	0.0238799428886510
%         counter=1;
%         Rveccounter=1;
%         while counter<=size(indices,1)
%             I=intersect(find(indices(:,1)==indices(counter,1)),find(indices(:,2)==indices(counter,2)));    %I=intersect(find(indices(:,end-(i-1)*2-1)==indices(counter,end-(i-1)*2-1)),find(indices(:,end-(i-1)*2)==indices(counter,end-(i-1)*2)))
%             counter=I(end)+1; 
%             R{i}(:,I)=temp(:,Rveccounter)*ones(1,length(I));            
%             Rveccounter=Rveccounter+1;
%         end
%     else
%         counter=1;
%         Rveccounter=1;                        %i=3时，end-(i-1)*2-1=end-1,倒数第二列；end-(i-1)*2=end,倒数第一列；i=2时，end-(i-1)*2-1=end-3，end-(i-1)*2=end-2针对的是indices矩阵前面的几列
%         while counter<=size(indices,1)        %indices(:,end-1)==indices(1,end-1)
%             I=intersect(find(indices(:,end-(i-1)*2-1)==indices(counter,end-(i-1)*2-1)),find(indices(:,end-(i-1)*2)==indices(counter,end-(i-1)*2)));%第一次循环，I={1，2}{3,4}{5,6};;第二次循环I={1}{2}{3}{4}{5}{6}
%             counter=I(end)+1;%counter = 3 
%             R{i}(:,I)=temp(:,Rveccounter)*ones(1,length(I));%将R(:,1)进行扩展，扩展成indices最后一列的一组。即level1的子节点个数，逐个扩展；存储在R{1}的各个列里面 ;;;       
%             Rveccounter=Rveccounter+1;                      %将R(:,2)进行扩展，将R(:,3)进行扩展
%         end                                                 %后续会存储到R{2}的各个列里面；将R{2}(:,1)=temp(:,1)*ones(1,1)
%     end
% end
% 
% 
% 
% % construct the core tensor
% S=zeros(prod(r),1);%以叶子节点的个数为行，构建一个列向量
% for k=1:length(nonzerosigmas)
%     temp=nonzerosigmas(k)*kron(R{d}(:,k),R{d-1}(:,k));%temp=nonzerosigmas(1)*kron(R{3}(:,1),R{2}(:,1)) nonzerosigmas(2)*kron(R{3}(:,2),R{2}(:,2))
%     for i=d-2:-1:1%这里d-2是因为已经计算了两个维度，或者说是
%         temp=kron(temp,R{i}(:,k));  %kron(temp,R{1}(:,1))  这里构造的是一个与原来的张量元素的个数相同的一个列向量kron(temp,R{2}(:,2))
%     end
%     S=S+temp;    
% end
% 
% 
% 
% tol=max(n)*eps(max(nonzerosigmas));%构造一个阈值，当S的元素小于这个阈值的时候，可以将这个元素默认为0
% S(abs(S)<tol)=0;
% S=reshape(S,r);
% 
% In multilinear algebra, there does not exist a general decomposition method for multi-way arrays (also known as N-arrays, higher-order arrays, or data-tensors) with all the properties of a matrix singular value decomposition (SVD). A matrix SVD simultaneously computes
% (a) a rank-R decomposition and
% (b) the orthonormal row/column matrices.
% These two properties can be captured separately by two different decompositions for multi-way arrays.
% Property (a) is extended to higher order by a class of closely related constructions known collectively as CP decomposition (named after the two most popular and general variants, CANDECOMP and PARAFAC). Such decompositions represent a tensor as the sum of the n-fold outer products of rank-1 tensors, where n is the dimension of the tensor indices.
% Property (b) is extended to higher order by a class of methods known variably as Tucker3, N-mode SVD, and N-mode principal component analysis (PCA). (This article will use the general term "Tucker decomposition".) These methods compute the orthonormal spaces associated with the different axes (or modes) of a tensor. The Tucker decomposition is also used in multilinear subspace learning as multilinear principal component analysis. This terminology was coined by P. Kroonenberg in the 1980s, but it was later called multilinear SVD and HOSVD (higher-order SVD) by L. De Lathauwer.
% Historically, much of the interest in higher-order SVDs was driven by the need to analyze empirical data, especially in psychometrics and chemometrics. As such, many of the methods have been independently invented several times, often with subtle variations, leading to a confusing literature. Abstract and general mathematical theorems are rare (though see Kruskal[1] with regard to the CP decomposition); instead, the methods are often designed for analyzing specific data types. The 2008 review article by Kolda and Bader[2] provides a compact summary of the history of these decompositions, and many references for further reading.
% The concept of HOSVD was carried over to functions by Baranyi and Yam via the TP model transformation [3] .[4] This extension led to the definition of the HOSVD based canonical form of tensor product functions and Linear Parameter Varying system models [5] and to convex hull manipulation based control optimization theory, see TP model transformation in control theories.