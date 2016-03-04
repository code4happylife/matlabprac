dt = 0.05;%在连续区间上取离散的计算点，设定采样间隔
t = 0:dt:5;%数值计算只能在有限的区间上，取有限个采样点
Ft = t.^2.*cos(t);%函数关系
Sx = dt*cumtrapz(Ft);%通过梯形近似计算积分
t(end -4:end)%打印出最后五个t的值
Sx(end -4:end)%打印出与最后五个t相对应的最后五个Sx的值
plot(t,Sx,'.k','MarkerSize',12)
xlabel('x'),ylabel('Sx'),grid on