function vec = matrix2vector(A)
sz=size(A);
num=sz(1)*sz(2);
vec=zeros(1,num);
count=1;
for j=1:sz(2)
    for i=1:sz(1)
        vec(count)=A(i,j);
        count=count+1;
    end
end
end

    