clear%清除所有的内存变量
a=12345.6789%数值标量
class(a)%对变量a的类别进行判断
a_s=size(a)%数值数组a的维度
b='S'%给变量b赋字符标量
class(b)%符号数组b的大小
b_s=size(b)
whos%观察变量a,b在内存中所占的字节

% numberandstring
% 
% a =
% 
%    1.2346e+04
% 
% 
% ans =
% 
% double
% 
% 
% a_s =
% 
%      1     1
% 
% 
% b =
% 
% S
% 
% 
% ans =
% 
% char
% 
% 
% b_s =
% 
%      1     1
% 
%   Name      Size            Bytes  Class     Attributes
% 
%   a         1x1                 8  double              
%   a_s       1x2                16  double              
%   ans       1x4                 8  char                
%   b         1x1                 2  char                
%   b_s       1x2                16  double       
