function tensor= vectorAndtensorIndexRelationship(vec,inputsize)
%��һ����������Ԫ��ȡ������ʹ�����չ����һ��������������̣����������Ĺ���
%����̽������������ÿһ���������һά����������Ĺ�ϵ���⡣
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

%size���������ĸ���ά�ȣ�����[2,3,4],����һ��size�����������ĳһ���ض���������
%����ÿ��ά�ȣ������Ӧ��index,Ҫ��������������ʾ��������������������size�����ʾ��
%size��һ��һ�У�sizedimensions�е�����������ȷ���ģ����ǣ�ÿһά�ȵ���ֵ��size(1),size(2),size(3),size(4),size(5)����ν���������ϣ���ȡÿһ��tensor����Ŀ
%tensor(1,2,3,4,5)=vec((5-1)*size(4)*size(3)*size(2)*size(1)+(4-1)*size(3)*size(2)*size(1)+...)
%�����������ڻ�
terms=[];

for n=1:sizedimension
    terms=[terms indicesFactorial(n)];%1*nȷ���Ե�
end
%����������Ǳ仯��
%��һ��������tensorindex{i}(:,j),i�����ֵ��sizedimension,j�����ֵ��size(1),j=1:size(1)������iһ��Ҫ���ϵ��������ƽ���j������ν
%�ڶ���������tensorindex{2}(:,k),i�����ֵ��sizedimension,k�����ֵ��size(2),k=1:size(2)
%������������tensorindex{3}(:,l),i�����ֵ��sizedimension,l�����ֵ��size(3),l=1:size(3)

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
