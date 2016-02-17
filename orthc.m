function O=orthc(A,varargin)
A(:,:,1)=[1 4 7 10;2 5 8 11;3 6 9 12];
A(:,:,2)=[ 13 16 19 22; 14 17 20 23; 15 18 21 24];%����һ��3������
n=size(A);%������ά��
r=zeros(1,length(n)-1);%zeros(1,2)
for i=1:length(n)-1
    r(i) = min(n(i),prod(n(i+1:end)));
end
totalsvd=1;
svdsperlevel=zeros(1,length(r));   % add length 1 for the first level
svdsperlevel(1)=1;  % first level
for i=2:length(r)
    svdsperlevel(i)=prod(r(1:i-1));
    totalsvd=totalsvd+svdsperlevel(i); %��ȷҪ�����svd����Щ��ΪʲôҪ������Щsvd������㷨�ĺ�����ʲô��
end
nleaf=prod(r);

U=cell(1,totalsvd);
S=cell(1,totalsvd);
V=cell(1,totalsvd);

[Ut St Vt]=svd(reshape(A,[n(1),prod(n(2:end))]));%��ԭ������������reshape�任��������չ���ɾ�����ö�Ӧ��svd
U{1}=Ut;
S{1}=diag(St);
V{1}=Vt;
counter=2; % this counter keeps track during the iterations which V{i} we are computing. This is a linear counter that counts breadth-first
whichvcounter=1;    % this counter keeps track during the iterations of which V we are taking the svd

for i=1:length(r)-1           % outer loop over the levels
    for j=1:prod(r(1:i))      % inner loop over the number of svds for this level 
        if rem(j,r(i)) == 0
            col=r(i);
        else
            col=rem(j,r(i));
        end
        [Ut St Vt]=svd(reshape(V{whichvcounter}(:,col),[n(i+1),prod(n(i+2:end))]));
        U{counter}=Ut;
        S{counter}=diag(St);
        V{counter}=Vt;
        counter=counter+1;
        if rem(j,length(S{whichvcounter}))==0
%             V{whichvcounter}=[];
            whichvcounter =  whichvcounter+1;
        end
    end
%     whichvcounter = whichvcounter+1;
end

Slevel=cell(1,length(r));   % cat each level singular values into 1 vector
counter=1;
for i=1:length(r),%����2��ѭ���ĳ��ֺ�Ƶ������Ҫ����ע�������
    for j=1:svdsperlevel(i),
        Slevel{i}=[Slevel{i}; S{counter}];
        counter=counter+1;
    end
end

for i=1:length(r),             % make all singular value vectors the same size (number of leaves)
    Slevel{i}=kron(Slevel{i}, ones(nleaf/length(Slevel{i}),1));
end

sigmas=ones(nleaf,1);         % output singular values at each leaf
for i=1:length(r),
    sigmas=sigmas.*Slevel{i};
end

if nargin==2
    tol=varargin{1};
else
    tol= max(n)*eps(sigmas(1));%��������һ����ֵ������������Щ��С�����ݣ�
    %�����ֵ������ĳ������ֵ�Ƿ������Ϊ0���򻯾���ı�ʾ��
    %The tensor train rank-1 decomposition exhibits
    %a singular value pro?le as with the SVD, 
    %allowing for a low-rank truncated series whenever the
    %singular value decay is prominent.
end

counter=1;
for i=1:svdsperlevel(end)
    %% for each svd of the last level
    
    % how many leaves per svd?
    leavespersvd=prod(r)/svdsperlevel(end);%����������tree�������ݽṹ
    indices=leave2ind((i-1)*leavespersvd+1,n);%���������leave2ind����������������������Ŀ����ʱ����ȷ    
    numberofUs=size(U{indices(1)},2);
    
    % check for numerically zero sigmas
    I=find(sigmas((i-1)*leavespersvd+1:i*leavespersvd)<tol);
    if ~isempty(I)
        for k=1:length(I)
            for j=1:length(r)-1
                O{j,counter}=U{indices( (length(r)-j)*2+1 )}(:,indices((length(r)-j)*2+2));
            end
            O{length(r),counter}=U{indices(1)}(:,I(k));
            O{length(r)+1,counter}=V{indices(1)}(:,I(k));
            counter =counter+1;
        end        
    end
    
    for k=1:numberofUs
        Vindices=1:n(end);
        if k <=n(end)
            Vindices(k)=[];
        end
        for l=1:length(Vindices)        % first we go down the tree up until the last level            
            for j=1:length(r)-1
                O{j,counter}=U{indices( (length(r)-j)*2+1 )}(:,indices((length(r)-j)*2+2));
            end
            O{length(r),counter}=U{indices(1)}(:,k);
            O{length(r)+1,counter}=V{indices(1)}(:,Vindices(l));
            counter =counter+1;
        end
    end    
end

end
% O=orthc(A) or O=orthc(A,tol)
% ----------------------------
% Computes the outer product vectors that span the orthogonal complement
% tensors.
%���㹹�������걸���������������
% O         =   cell, each column of O contains outer product vectors that form
%               a tensor orthogonal to A,
%
% A         =   array, d-way array,
%
% tol       =   scalar, optional tolerance to decide whether computed singular
%               values are numerically zero. Default value = max(size(A))*eps(sigmas(1)).
%