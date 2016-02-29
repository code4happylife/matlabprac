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

d=length(n);                            % order of the tensor��ÿһ�γ������������Բ�ͬ�Ľ������������з���
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
endI=endI(1:end-2); %endI�ĳ��ȱ�ʵ�ʵ�ά����2����Ҫ����������ÿһ�������svd����

% k(i) containst the offset of the ith term with respect to endI for
% level d-i, initalize to sigmaI
k=sigmaI;
% sigmaI��¼�������з����Ҷ�ӽڵ�����Ӧ������
% l tells us which element of S{k(i)}(:,l(i)) we need to choose
l=zeros(1,length(sigmaI));%k,l��S�Ĺ�ϵ����

indices=zeros(length(sigmaI),2*(length(n)-1));%��С��length(n)-1��ά�Ƚ��з���

for j=1:length(sigmaI)%��ÿһ��Ҷ�ӽڵ㿪ʼ����������ÿһ��Ҷ�ӽڵ��������ʼ�ɻ�
    
    % determine offset for each node
    if mod(k(j),r(end))==0%r(end)=2,������������ƫ�������з���
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
        % determine offset(λ��) for each node 
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
%�����indices����Ľ���Ǽ�����
% ֻ�ǣ���Ȼ��̫������ľ��庬�壬��Ҫ��һ���Ķ�svd�ֽ��qr�ֽ����ѧϰ�о���������Ҫ�Ժ����Ǿ�����������indices�и�����İѿ�
% indices =
% 
%      2     1     1     1
%      2     2     1     1
%      3     1     1     2
%      3     2     1     2
%      4     1     1     3
%      4     2     1     3
%���������������ж���indices�ĸ����¼����ÿ������ֵ�ڲ�ͬ��S�����е�����ֵ��������Щ����ֵ���Ը��õ�access��Щ����ֵ���������ڶ��ڹ�С������ֵ������