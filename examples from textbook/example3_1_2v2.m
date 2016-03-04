dt = 0.05;%�趨�������
t = 0:dt:4;
Ft = exp(-sin(t));%����ȷ���ĺ�����ϵ��������ÿ�����������Ӧ�ĺ���ֵ
Sx = dt*cumtrapz(Ft);
Sx(end)
plot(t,Ft,'*r','MarkerSize',4)
hold on
plot(t,Sx,'.k','MarkerSize',15)
hold off
xlabel('x')
legend('Ft','Sx')