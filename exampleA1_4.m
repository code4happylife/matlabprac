clear%清除内存中的所有变量
a=2;%设置衰减系数
w=3;%设置振荡频率
t=0:0.01:10;%取自变量采样数组
y=exp(-a*t).*sin(w*t);%计算函数值，产生函数数组
[y_max,i_max]=max(y);%找到最大元素的位置
t_text=['t=',num2str(t(i_max))];%生成最大值点的横坐标字符串
y_text=['y=',num2str(y_max)];%生成最大值点的纵坐标字符串
max_text=char('maximum',t_text,y_text);%生成标志最大值点的三行字符串
tit=['y=exp(-',num2str(a),'t)*sin(',num2str(w),'t)'];%生成标志图名用的字符串
plot(t,zeros(size(t)),'k')%画纵坐标为0的基准线
hold on
plot(t,y,'b')%用蓝色画y(t)曲线
plot(t(i_max),y_max,'r.','MarkerSize',20)%用大红色点标最大值点
text(t(i_max)+0.3,y_max+0.05,max_text)%在图上书写最大值点的数据值
title(tit),xlabel('t'),ylabel('y')%书写图名，横坐标名，纵坐标名
hold off
