clear%����ڴ��е����б���
a=2;%����˥��ϵ��
w=3;%������Ƶ��
t=0:0.01:10;%ȡ�Ա�����������
y=exp(-a*t).*sin(w*t);%���㺯��ֵ��������������
[y_max,i_max]=max(y);%�ҵ����Ԫ�ص�λ��
t_text=['t=',num2str(t(i_max))];%�������ֵ��ĺ������ַ���
y_text=['y=',num2str(y_max)];%�������ֵ����������ַ���
max_text=char('maximum',t_text,y_text);%���ɱ�־���ֵ��������ַ���
tit=['y=exp(-',num2str(a),'t)*sin(',num2str(w),'t)'];%���ɱ�־ͼ���õ��ַ���
plot(t,zeros(size(t)),'k')%��������Ϊ0�Ļ�׼��
hold on
plot(t,y,'b')%����ɫ��y(t)����
plot(t(i_max),y_max,'r.','MarkerSize',20)%�ô��ɫ������ֵ��
text(t(i_max)+0.3,y_max+0.05,max_text)%��ͼ����д���ֵ�������ֵ
title(tit),xlabel('t'),ylabel('y')%��дͼ����������������������
hold off
