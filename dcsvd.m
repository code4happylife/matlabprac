function [U,S,V,sigmas]=dcsvd(A)
% [U,S,V,sigmas]=dcsvd(A)
% -------------------------
% Decomposes an arbitrary tensor A into
% a linear combination of orthonormal rank-1 terms. Returns the orthogonal
% vectors U,V and singular values S from each of the SVDs in the tree.
% This version of tensor SVD determines numerically zero intermediate singular
% values and does not compute those branches of the TTr1-tree. This might
% result in a huge reduction of run time when the number of nonzero terms is small.
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
A(:,:,1)=[1 4 7 10;2 5 8 11;3 6 9 12];
A(:,:,2)=[ 13 16 19 22; 14 17 20 23; 15 18 21 24];
tensorSize=size(A);
r=zeros(1,length(tensorSize)-1);
for i=1:length(tensorSize)-1
    r(i) = min(tensorSize(i),prod(tensorSize(i+1:end)));
end

totalsvd=1;
svdsperlevel=zeros(1,length(r));   % add length 1 for the first level
svdsperlevel(1)=1;  % first level
for i=2:length(r)
    svdsperlevel(i)=prod(r(1:i-1));
    totalsvd=totalsvd+svdsperlevel(i); 
end
nleaf=prod(r);

U=cell(1,totalsvd);
S=cell(1,totalsvd);
V=cell(1,totalsvd);

[Ut St Vt]=svd(reshape(A,[tensorSize(1),prod(tensorSize(2:end))]),'econ');
U{1}=Ut;
S{1}=diag(St);
V{1}=Vt;
counter=2; % this counter keeps track during the iterations which V{i} we are computing. This is a linear counter that counts breadth-first
whichvcounter=1;    % this counter keeps track during the iterations of which V we are taking the svd
sigmas=kron(S{1}, ones(nleaf/length(S{1}),1));

tol=prod(tensorSize)*eps(max(S{whichvcounter}));

for i=1:length(r)-1           % outer loop over the levels
    Slevel=[];
    for j=1:prod(r(1:i))      % inner loop over the number of eigs for this level 
        if rem(j,r(i)) == 0
            col=r(i);
        else
            col=rem(j,r(i));
        end
        if ~isempty(S{whichvcounter}) && S{whichvcounter}(col) > tol
			[Ut St Vt]=svd(reshape(V{whichvcounter}(:,col),[tensorSize(i+1),prod(tensorSize(i+2:end))]),'econ');
			U{counter}=Ut;
			S{counter}=diag(St);
			V{counter}=Vt;
			Slevel=[Slevel;S{counter}];            
        else
            Slevel=[Slevel;zeros(min([tensorSize(i+1),prod(tensorSize(i+2:end))]),1)];
        end
        counter=counter+1;
        if rem(j,tensorSize(i))==0
            whichvcounter =  whichvcounter+1;
        end
    end
    Slevel=kron(Slevel, ones(nleaf/length(Slevel),1));
    sigmas=sigmas.*Slevel;
end

end
