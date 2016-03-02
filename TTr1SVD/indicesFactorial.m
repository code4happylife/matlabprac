function newindex =indicesFactorial(index)
n=size(tensor);
newindex=1;
if index~=1
    for i=index-1:-1:1
        newindex=newindex*n(i);
    end
else
    newindex=1;
end
