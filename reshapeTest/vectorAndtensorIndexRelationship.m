function tensor= vectorAndtensorIndexRelationship(vec,inputsize)
%将一个行向量的元素取出来，使它最终构造出一个张量，这个过程，是张量化的过程
%这里探究的是张量的每一项的索引与一维数组的索引的关系问题。
% vec=[1 2 3 4 5 6 7 8 9 10 11 12];
% size=[2 3 2];
global tensorsize;
tensorsize = inputsize;
sizedimension=length(tensorsize);
tensor=zeros(tensorsize);
% tensor=vec(1:sizedimension)

% tensorindex=cell(1,sizedimension);
% for i=1:sizedimension
%     tensorindex{i}=(1:inputsize(i));
% end

%size代表张量的各个维度，比如[2,3,4],依据一个size向量怎样表达某一个特定的张量？
%对于每个维度，计算对应的index,要将张量的索引表示出来，张量的索引是由size数组表示的
%size是一个一行，sizedimensions列的向量，这是确定的，但是，每一维度的数值，size(1),size(2),size(3),size(4),size(5)则如何进行排列组合，存取每一个tensor的项目
%tensor(1,2,3,4,5)=vec((5-1)*size(4)*size(3)*size(2)*size(1)+(4-1)*size(3)*size(2)*size(1)+...)
%索引看成是内积
terms=[];

for n=1:sizedimension
    terms=[terms indicesFactorial(n)];%1*n确定性的
end
%这里的索引是变化的
%第一个索引是tensorindex{i}(:,j),i的最大值是sizedimension,j的最大值是size(1),j=1:size(1)，其中i一定要不断的往后面推进，j则无所谓
%第二个索引是tensorindex{2}(:,k),i的最大值是sizedimension,k的最大值是size(2),k=1:size(2)
%第三个索引是tensorindex{3}(:,l),i的最大值是sizedimension,l的最大值是size(3),l=1:size(3)

% count =0;
% cell tensor_index = cell(1,sizedimension);
% for iter = 1:sizedimension
%    tensor_index=[tensor_index tensorindex{iter}];
% end
tensorindex.indices=cell(1,sizedimension);
tensorindex.actual=[];
for i=1:sizedimension
    
    matrix=1:inputsize(i)
    col=cell(1,numel(matrix));
    tensorindex.indices{i}=matrix;  
    
    for j=1:col
        tensorindex.actual=[tensorindex.actual col]        
    end
    
end
tensor(tensorindex(i))

tensorindex()
for i=1:sizedimension
    j=1:size(inputsize(i));
    cnt=cnt+tensorindex(i).iteration(j)*terms(i);
    tensor(tensorindex{:})=vec(cnt);
end


    function newindex =indicesFactorial(index)
        newindex=1;
        if index~=1
            for i=index-1:-1:1
                newindex=newindex*tensorsize(i);
            end
        else
            newindex=1;
        end
    end
end
