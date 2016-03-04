function matrix = vector2matrix(vec,size)
matrix=zeros(size);
count=1;
for j=1:size(2)
    for i=1:size(1)
        matrix(i,j)=vec(count);
        count=count+1;
    end
end
end