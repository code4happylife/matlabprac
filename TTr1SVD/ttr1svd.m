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


%n��¼������ά�ȣ�����3*4*5*6*7������һ������
n=size(A);
r=zeros(1,length(n)-1);%rҲ��һ�����飬ֻ������ά�Ȼ����һ������ȫ����ʼ��Ϊ0.
for i=1:length(n)-1%ͬ�������������ĳ��ȣ��ֱ�����ǽ��и�ֵ��
    r(i) = min(n(i),prod(n(i+1:end)));%��ֵ��ʱ��ע�⣬r��ֵ�ֱ����n�ĸ���ά��ֵ�У���С��length��n��-1����r(i)�������ǣ��Ӷ��㿪ʼ��ÿһ��Ľڵ��������ӽڵ����Ŀ��
end%r(i)���Ǻ�������ģ��ֱ�ָ���˽�С��n-1��ά��ֵ


totalsvd=1;%������Ҫ�����svd����Ŀ
svdsperlevel=zeros(1,length(r));   % add length 1 for the first level����length(r):����r�ĳ���ָ�������νṹ�Ĳ�������3����������r�ĳ���Ϊ2����Ӧ���������ݽṹ����2������������Ҷ�ӽڵ��
svdsperlevel(1)=1;  % first level 
for i=2:length(r)
    svdsperlevel(i)=prod(r(1:i-1));%svdsperlevel����Ԫ�صĸ�ֵ������r�ĸ�����svdsperlevel�ĸ�����������ÿһ���ж��ٸ��ڵ���Ҫ���㣻��Ӧ��ȥ����ѭ���ؼ���svd
    totalsvd=totalsvd+svdsperlevel(i);%totalsvd��svdsperlevel(2)��ʼ���������� 
end
nleaf=prod(r);%Ҷ�ӽڵ����������r��ÿ��Ԫ�ؽ��г˻����㣻�����������һ��Ľڵ����

U=cell(1,totalsvd);%��Ԫ���ݣ�totalsvd�����е�svd����Ŀ
S=cell(1,totalsvd);
V=cell(1,totalsvd);

[Ut St Vt]=svd(reshape(A,[n(1),prod(n(2:end))]),'econ');%��һ�μ���svd�����������а��н���reshape����n(1)Ϊ�У���prod(n(2:end))Ϊ�У�������������svd����
U{1}=Ut;%��Ut�Ľ���洢��U{1}��Ut��n(1)*n(1)
S{1}=diag(St);%��St�ĶԽ���洢��S{1}��St��n(1)*prod(n(2:end))
V{1}=Vt;%Vt��prod(n(2:end))*prod(n(2:end))
counter=2; % this counter keeps track during the iterations which V{i} we are computing. This is a linear counter that counts breadth-first
whichvcounter=1;    % this counter keeps track during the iterations of which V we are taking the svd
sigmas=kron(S{1}, ones(nleaf/length(S{1}),1));%length(S{1})=n(1)��nleaf=prod(r)��r(i) = min(n(i),prod(n(i+1:end)))

for i=1:length(r)-1           % outer loop over the levels����ǰ������νṹ�Ļ����ϣ��ٱ���һ�㣻length(r)-1=length(n)-2��3�������������һ�㣻4�������������2�㣻���һ���Ҷ�ӽڵ㲻��Ҫ����
    Slevel=[];
    for j=1:prod(r(1:i))      % inner loop over the number of svds for this level��svdsperlevel(i)=prod(r(1:i-1))��ʵ�������prod��r(1:i)������ָsvdsperlevel(2)����ÿһ���ж��ٸ��ڵ���Ҫ���㣻j=[1:nodes(2)]��
        if rem(j,r(i)) == 0   %r(i)�������ǣ��Ӷ��㿪ʼ��ÿһ��Ľڵ��������ӽڵ����Ŀ
            col=r(i);         
        else
            col=rem(j,r(i)); %ƫ�����ļ��㣻
        end
        [Ut St Vt]=svd(reshape(V{whichvcounter}(:,col),[n(i+1),prod(n(i+2:end))]),'econ');
        %V{1}��������prod(n(2:end))������Ϊ8��r(1)=3����������ֵ�ֽ���ص㣬���ֻ��Ҫ��V��ǰm=n(1)=3=r(1)�����������к������㼴�ɣ�col
        %������r(i)
        U{counter}=Ut;%Ut��ά��=[n(i+1)*n(i+1)]
        S{counter}=diag(St);%St��ά��=[n(i+1),prod(n(i+2):end)]��diag(St)=n(i+1)��Ԫ��
        V{counter}=Vt;%Vt��ά��=[prod(n(i+2):end)*prod(n(i+2):end)]
        Slevel=[Slevel;S{counter}];
        counter=counter+1;
        if rem(j,length(S{whichvcounter}))==0%length(S{1})=n(1)
            whichvcounter =  whichvcounter+1;
        end
    end
    Slevel=kron(Slevel, ones(nleaf/length(Slevel),1));%�����nleaf�����һ��Ľڵ���Ŀ������Ŀ���ǰ�һ�ж���չ��Ҷ�ӽڵ���ô�������飻
    
    
    sigmas=sigmas.*Slevel;%��һ����õ�����ֵ���������һ��չ����õ�����ֵ�����ۼƳ˻�����
end
end
