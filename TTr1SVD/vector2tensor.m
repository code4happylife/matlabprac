function tensor= vector2tensor(vec,size)






sizedimension=length(size);
tensor=zeros(size);

tensorindex=cell(1,sizedimension);
for i=1:sizedimension
    tensorindex{i}=1:size(i);
end

for i=1:sizedimension
   
tensor(tensorindex{i}(:j))=vec()

count=1;

n=size;%n是一个数组
% for 
% count=(sn-1)*size(end-1)*size(end-2)*size(1)+(sn-1-1)(dn-2)...(d1)+...+(s2-1)(d1)+s1
% tensor()=vec(count)



tensor(size(n))=zeros(tensor(size(n)));
arr1=[];
for i=sizedimension:-1:1
    newindex =indicesFactorial(i)
    arr1=[arr1 newindex];   
end

iter=length(n);

arr2=[];
arr2=[]

arr1[].*arr2[]
tensor[]=vec(count);

for 

    for k=1:size(3)
 
                
                count=(k-1)*size(2)*size(1)+(j-1)*size(1)+i;
                tensor(i,j,k)=vec(count);
                
            end
        end
    end
else 
    for l=size(4)
        for k=1:size(3)
            for j=1:size(2)
                for i=1:size(1)
                    tensor(i,j,k,l)=vec(count);
                    count=count+1;
                end
            end
        end
    end
end


