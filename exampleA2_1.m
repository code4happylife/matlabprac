clear
C_str='这是胞元数组创建算例1';%产生字符串
R=reshape(1:9,3,3);%产生（3*3）实数矩阵R
Cn=[1+2i];%产生复数标量
S_sym=sym('sin(-3*t)*exp(-t)');%产生符号函数量
%创建胞元数组方法之一
B{1,1}=C_str;
B{1,2}=R;
B{2,1}=Cn;
B{2,2}=S_sym;
%胞元的援引
a=B(1,2)%注意这里使用圆括号
class(a)%a是胞元
%胞元内容的援引
b=B{1,2}%注意这里使用花括号
class(b)%b是双精度矩阵

% exampleA2_1
% 
% a = 
% 
%     [3x3 double]
% 
% 
% ans =
% 
% cell
% 
% 
% b =
% 
%      1     4     7
%      2     5     8
%      3     6     9
% 
% 
% ans =
% 
% double
% B{2,2}
%  
% ans =
%  
% -sin(3*t)*exp(-t)