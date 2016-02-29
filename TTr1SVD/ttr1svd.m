function [U,S,V,sigmas]=ttr1svd(A)
% [U,S,V,sigmas]=ttr1svd(A)
% -------------------------
% Tensor Train rank-1 singular value decomposition. Decomposes an arbitrary tensor A into
% a linear combination of orthonormal rank-1 terms. Returns the orthogonal
% vectors U,V and singular values S from each of the SVDs in the TTr1 tree.
% Use the function getAtilde.m to obtain a low rank approximation using the
% U,S,V obtained from this function.
%
% U         =   cell, contains the U vectors of each of the SVDs in the
%               tree,
%
% S         =   cell, contains the singular values S of each of the SVDs
%               in the tree,
%
% V         =   cell, contains the V vectors of each of the SVDs in the
%               tree,
%
% sigmas    =   vector, contains the final singular values in the linear
%               combination of rank-1 terms.


%n记录张量的维度，比如3*4*5*6*7，它是一个数组
n=size(A);
r=zeros(1,length(n)-1);%r也是一个数组，只是它的维度会减少一个，且全部初始化为0.
for i=1:length(n)-1%同理，根据这个数组的长度，分别对它们进行赋值。
    r(i) = min(n(i),prod(n(i+1:end)));%赋值的时候注意，r的值分别等于n的各个维度值中，较小的length（n）-1个；r(i)的意义是，从顶点开始，每一层的节点生出的子节点的数目；
end%r(i)都是很有意义的，分别指代了较小的n-1个维度值


totalsvd=1;%整个需要计算的svd的数目
svdsperlevel=zeros(1,length(r));   % add length 1 for the first level这里length(r):数组r的长度指的是树形结构的层数，即3阶张量，则r的长度为2；对应的树型数据结构，有2层树；不包括叶子节点层
svdsperlevel(1)=1;  % first level 
for i=2:length(r)
    svdsperlevel(i)=prod(r(1:i-1));%svdsperlevel各个元素的赋值；根据r的各项，求得svdsperlevel的个数，即，树每一层有多少个节点需要计算；对应的去迭代循环地计算svd
    totalsvd=totalsvd+svdsperlevel(i);%totalsvd从svdsperlevel(2)开始计数，叠加 
end
nleaf=prod(r);%叶子节点的总数，将r的每个元素进行乘积运算；这是树的最后一层的节点个数

U=cell(1,totalsvd);%胞元数据，totalsvd是所有的svd的数目
S=cell(1,totalsvd);
V=cell(1,totalsvd);

[Ut St Vt]=svd(reshape(A,[n(1),prod(n(2:end))]),'econ');%第一次计算svd，将张量进行按列进行reshape，以n(1)为行，以prod(n(2:end))为列；对这个矩阵进行svd计算
U{1}=Ut;%将Ut的结果存储到U{1}；Ut：n(1)*n(1)
S{1}=diag(St);%将St的对角项存储在S{1}；St：n(1)*prod(n(2:end))
V{1}=Vt;%Vt：prod(n(2:end))*prod(n(2:end))
counter=2; % this counter keeps track during the iterations which V{i} we are computing. This is a linear counter that counts breadth-first
whichvcounter=1;    % this counter keeps track during the iterations of which V we are taking the svd
sigmas=kron(S{1}, ones(nleaf/length(S{1}),1));%length(S{1})=n(1)；nleaf=prod(r)；r(i) = min(n(i),prod(n(i+1:end)))

for i=1:length(r)-1           % outer loop over the levels，在前面的树形结构的基础上，少遍历一层；length(r)-1=length(n)-2；3阶张量，则遍历一层；4阶张量，则遍历2层；最后一层的叶子节点不需要遍历
    Slevel=[];
    for j=1:prod(r(1:i))      % inner loop over the number of svds for this level；svdsperlevel(i)=prod(r(1:i-1))；实际这里的prod（r(1:i)）就是指svdsperlevel(2)；树每一层有多少个节点需要计算；j=[1:nodes(2)]；
        if rem(j,r(i)) == 0   %r(i)的意义是，从顶点开始，每一层的节点生出的子节点的数目
            col=r(i);         
        else
            col=rem(j,r(i)); %偏移量的计算；
        end
        [Ut St Vt]=svd(reshape(V{whichvcounter}(:,col),[n(i+1),prod(n(i+2:end))]),'econ');
        %V{1}的列数：prod(n(2:end))；假设为8；r(1)=3；由于奇异值分解的特点，因此只需要对V的前m=n(1)=3=r(1)个列向量进行后续计算即可；col
        %最大等于r(i)
        U{counter}=Ut;%Ut的维度=[n(i+1)*n(i+1)]
        S{counter}=diag(St);%St的维度=[n(i+1),prod(n(i+2):end)]；diag(St)=n(i+1)个元素
        V{counter}=Vt;%Vt的维度=[prod(n(i+2):end)*prod(n(i+2):end)]
        Slevel=[Slevel;S{counter}];
        counter=counter+1;
        if rem(j,length(S{whichvcounter}))==0%length(S{1})=n(1)
            whichvcounter =  whichvcounter+1;
        end
    end
    Slevel=kron(Slevel, ones(nleaf/length(Slevel),1));%这里的nleaf是最后一层的节点数目；最后的目的是把一切都扩展成叶子节点那么长的数组；
    
    
    sigmas=sigmas.*Slevel;%第一次求得的奇异值，与后续进一步展开求得的奇异值进行累计乘积运算
end
end
