dt = 0.05;%设定采样间隔
t = 0:dt:4;
Ft = exp(-sin(t));%按照确定的函数关系，计算与每个采样点相对应的函数值
Sx = dt*cumtrapz(Ft);
Sx(end)
plot(t,Ft,'*r','MarkerSize',4)
hold on
plot(t,Sx,'.k','MarkerSize',15)
hold off
xlabel('x')
legend('Ft','Sx')