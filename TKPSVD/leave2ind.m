function indices=leave2ind(sigmaI,n)
% indices=leave2ind(sigmaI,n)
% ---------------------------
% Returns the indices of the cells and column vectors of the U,V
% cells/vectors that are required to reconstruct the rank-1 terms
% determined by sigmaI.返回胞元U和V的列向量的索引.
%
% indices   =   matrix, each row contains pairs of cell/column indices
%
% sigmaI    =   vector, contains indices of rank-1 terms that need to be
%               used in the reconstruction,
%
% n         =   vector, dimensions of the original tensor A.


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

%last endI points to last SVD of previous level,这里的endI的变换至关重要       
endI=endI(1:end-2); 

% k(i) containst the offset of the ith term with respect to endI for
% level d-i, initalize to sigmaI
k=sigmaI;
% l tells us which element of S{k(i)}(:,l(i)) we need to choose
l=zeros(1,length(sigmaI));
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

for i=1:length(endI)%经过之前的变换，endI的长度已经变为0了，所以endI究竟想做什么呢？    
    for j=1:length(sigmaI)
        % determine offset for each node
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
    % determine offset for each node
    if mod(k(j),r(end))==0
        l(j)=r(end);
    else
        l(j)=mod(k(j),r(end));
    end    
    indices(j,Tcounter)=1;
    indices(j,Tcounter+1)=l(j); 
end

end