a0 = [0.2,pi/2,-2,sin(pi/5),-exp(-3)]
a1 = 1:6
a2 = 0:pi/4:pi
a3 = 1:-0.1:0
b1 = linspace(0,pi,4)%相当于0:pi/3:pi
b2 = logspace(0,3,4)%创建数组[10^0 10^1 10^2 10^3]
c1 = [2 pi/2 sqrt(3) 3+5i]
rng default
c2=rand(1,5)%产生(1*5)的均匀分布随机数组
x1=(1:6)'
x2=linspace(0,pi,4)'
y1=rand(5,1)
z1=[2;pi/2;sqrt(3);3+5i]