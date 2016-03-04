function tensor=vecTotensor(vec,inputsize)
% inputsize=[2,3,2,2];
% vec=zeros(1,prod(inputsize))
tensor=zeros(inputsize);
% for i=1:prod(inputsize)
% vec(i)=i;
% end
tensor(:)=vec(1:numel(vec));
end

%According to the mentor's advice,just use colon noation in Matlab can
%easily solve my problem of tranform a vector into a tensor with the same
%functionality as what the reshape command does in Matlab.




