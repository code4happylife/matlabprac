function [Q,R]=qrgv(A)
% ����Givens�任��������A�ֽ�ΪA=QR������QΪ��������RΪ��������
%
% ����˵��
% A����Ҫ����QR�ֽ�ķ���
% Q���ֽ�õ�����������
% R���ֽ�õ�����������
%
% ʵ��˵��
A=[-12 3 3;3 1 -2;3 -2 7];
% [Q,R]=qr(A) % ����MATLAB�Դ���QR�ֽ⺯��������֤
% [q,r]=qrgv(A) % ���ñ���������QR�ֽ�
% q*r-A % ��֤ A=QR
% q'*q % ��֤q��������
% norm(q) % ��֤q�ı�׼����������������1
%
% ���Դ�������֪ʶ
% 1.B=P*A*inv(P)����A��B���ƣ����ƾ��������ͬ������ֵ
% 2.Q*Q'=I����QΪ����������������ĳ˻���Ϊ��������
n=size(A,1);
R=A;
Q=eye(n);
for i=1:n-1
    for j=2:n-i+1
        x=R(i:n,i);
        rt=givens(x,1,j);%In numerical linear algebra, a Givens rotation is a rotation in the plane spanned by two coordinates axes. 
        %Givens rotations are named after Wallace Givens, 
        %who introduced them to numerical analysts in the 1950s while he was working at Argonne National Laboratory.
        r=blkdiag(eye(i-1),rt);%Construct block diagonal matrix from input arguments
        
        Q=Q*r';
        R=r*R;
    end
end
