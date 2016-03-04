% n-2个正交相似变换的序列
% 本模块演示通过正交相似变换A = Q H QT将矩阵A规约为上Hessenberg型矩阵.
% 矩阵H是上Hessenberg矩阵，即它的第一条次对角线下方的元素都是零，Q是正交矩阵.
% 变成上Hessenberg矩阵的规约过程是使用QR迭代计算矩阵特征值和特征向量的前提步骤，它通过运用Hoseholder变换逐列地消去矩阵的部分元素来实现.
% 正交相似变换保持矩阵的对称性，所以如果初始矩阵A是 对称的，那么结果矩阵将是三对角矩阵，这种情况的结果矩阵常用T, 而不是H来表示.
A=[1 2 3;4 5 6;7 8 9];
n=size(A,1);
for k = 1:n-2
u = A(:,k);
u(1:k) = 0;
sigma = norm(u);
if sigma ~= 0
    
if u(k+1) < 0, sigma = -sigma; 
end
u(k+1) = u(k+1) + sigma;
rho = 1/(sigma*u(k+1));
v = rho*A*u;
w = rho*(u'*A)';
gamma = rho/2*u'*v;
v = v - gamma*u;
w = w - gamma*u;
A = A - v*u' - u*w';
A(k+2:n,k) = 0;
end

end