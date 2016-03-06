function [indices,r,endI]=testleave2ind(sigmaI,n)
d=length(n);                            % order of the tensor
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
k=sigmaI;%sigmaI记录的是所有非零的叶子节点所对应的索引
% l tells us which element of S{k(i)}(:,l(i)) we need to choose
l=zeros(1,length(sigmaI));%k,l和S的关系密切

indices=zeros(length(sigmaI),2*(length(n)-1));

for j=1:length(sigmaI)
    
    % determine offset for each node
    if mod(k(j),r(end))==0
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