function indices=leave2ind(sigmaI,n)
% indices=leave2ind(sigmaI,n)
% ---------------------------
% Returns the indices of the cells and column vectors of the U,V
% cells/vectors that are required to reconstruct the rank-1 terms
% determined by sigmaI.
%
% indices   =   matrix, each row contains pairs of cell/column indices
%
% sigmaI    =   vector, contains indices of rank-1 terms that need to be
%               used in the reconstruction,
%
% n         =   vector, dimensions of the original tensor A.

% Reference
% ---------
%
% A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms
% http://arxiv.org/abs/1407.1593
%
% 2014, Kim Batselier, Haotian Liu, Ngai Wong

d=length(n);                            % order of the tensor，每一段程序都有这个，针对不同的阶数的张量进行分析
if d==2
    % matrix case
    indices=ones(length(sigmaI),2);
    indices(:,2)=[1:length(sigmaI)]';
    return
end

r=ones(1,d);                            % nodes are paired in groups of r(i) on level i-1
nodesperlevel=ones(1,d);                % level i-1 contains nodesperlevel(i) nodes
endI=ones(1,d);                         % endI(i) is the index of the last SVD of level i-1
for i=2:d
    r(i)=min(n(i-1),prod(n(i:end)));
    nodesperlevel(i)=prod(r(1:i));
    endI(i)=sum(nodesperlevel(1:i));
end

%last endI points to last SVD of previous level       
endI=endI(1:end-2); %endI的长度比实际的维度少2；需要遍历的树的每一层的最大的svd索引

% k(i) containst the offset of the ith term with respect to endI for
% level d-i, initalize to sigmaI
k=sigmaI;
% sigmaI记录的是所有非零的叶子节点所对应的索引
% l tells us which element of S{k(i)}(:,l(i)) we need to choose
l=zeros(1,length(sigmaI));%k,l和S的关系密切

indices=zeros(length(sigmaI),2*(length(n)-1));%较小的length(n)-1个维度进行分析

for j=1:length(sigmaI)%从每一个叶子节点开始分析，根据每一个叶子节点的索引开始干活
    
    % determine offset for each node
    if mod(k(j),r(end))==0%r(end)=2,结合它们自身的偏移量进行分析
        k(j)=k(j)/r(end);
        l(j)=r(end);
    else
        l(j)=mod(k(j),r(end));
        k(j)=ceil(k(j)/r(end));
    end
    
    indices(j,1)=endI(end)+k(j);
    indices(j,2)=l(j);    
end

endI=endI(1:end-1);
r=r(1:end-1);
Tcounter=3;
for i=1:length(endI)
    
    for j=1:length(sigmaI)
        % determine offset(位移) for each node 
        if mod(k(j),r(end))==0
            k(j)=k(j)/r(end);
            l(j)=r(end);
        else
            l(j)=mod(k(j),r(end));
            k(j)=ceil(k(j)/r(end));
        end
        
        indices(j,Tcounter)=endI(end)+k(j);
        indices(j,Tcounter+1)=l(j); 
                
    end    
    Tcounter=Tcounter+2;
    endI=endI(1:end-1);
    r=r(1:end-1);
end

for j=1:length(sigmaI)
    
    %% determine offset for each node
    if mod(k(j),r(end))==0
        l(j)=r(end);
    else
        l(j)=mod(k(j),r(end));
    end
    
    indices(j,Tcounter)=1;
    indices(j,Tcounter+1)=l(j); 
end

end
% [2,1,1,1;2,2,1,1;3,1,1,2;3,2,1,2;4,1,1,3;4,2,1,3]
%计算的indices矩阵的结果是简明的
% 只是，依然不太清楚它的具体含义，需要进一步的对svd分解和qr分解进行学习研究，并且需要对后续是具体怎样利用indices有个整体的把控
% indices =
% 
%      2     1     1     1
%      2     2     1     1
%      3     1     1     2
%      3     2     1     2
%      4     1     1     3
%      4     2     1     3
%经过分析，初步判定，indices的各项记录的是每个奇异值在不同的S矩阵中的索引值。根据这些索引值可以更好的access这些奇异值，包括后期对于过小的奇异值的舍弃