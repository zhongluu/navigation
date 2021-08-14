function [ MagDataOut ] = MagDataFitting( mag )
%% ���������������ϣ���С���˷���
%����ΰ
%�ο���ַ��https://blog.csdn.net/HJ199404182515/article/details/59480954?utm_source=app
%%
%mag:num_points*3
%%
x=mag(:,1);
y=mag(:,2);
z=mag(:,3);
%�ռ������������㷨
num_points = length(x);
%һ����ͳ��ƽ��
x_avr = sum(x)/num_points;
y_avr = sum(y)/num_points;
z_avr = sum(z)/num_points;
%������ͳ��ƽ��
xx_avr = sum(x.*x)/num_points;
yy_avr = sum(y.*y)/num_points;
zz_avr = sum(z.*z)/num_points;
xy_avr = sum(x.*y)/num_points;
xz_avr = sum(x.*z)/num_points;
yz_avr = sum(y.*z)/num_points;
%������ͳ��ƽ��
xxx_avr = sum(x.*x.*x)/num_points;
xxy_avr = sum(x.*x.*y)/num_points;
xxz_avr = sum(x.*x.*z)/num_points;
xyy_avr = sum(x.*y.*y)/num_points;
xzz_avr = sum(x.*z.*z)/num_points;
yyy_avr = sum(y.*y.*y)/num_points;
yyz_avr = sum(y.*y.*z)/num_points;
yzz_avr = sum(y.*z.*z)/num_points;
zzz_avr = sum(z.*z.*z)/num_points;
%�Ĵ���ͳ��ƽ��
yyyy_avr = sum(y.*y.*y.*y)/num_points;
zzzz_avr = sum(z.*z.*z.*z)/num_points;
xxyy_avr = sum(x.*x.*y.*y)/num_points;
xxzz_avr = sum(x.*x.*z.*z)/num_points;
yyzz_avr = sum(y.*y.*z.*z)/num_points;

%����������Է��̵�ϵ������
A0 = [yyyy_avr yyzz_avr xyy_avr yyy_avr yyz_avr yy_avr;
     yyzz_avr zzzz_avr xzz_avr yzz_avr zzz_avr zz_avr;
     xyy_avr  xzz_avr  xx_avr  xy_avr  xz_avr  x_avr;
     yyy_avr  yzz_avr  xy_avr  yy_avr  yz_avr  y_avr;
     yyz_avr  zzz_avr  xz_avr  yz_avr  zz_avr  z_avr;
     yy_avr   zz_avr   x_avr   y_avr   z_avr   1;];
%����������
b = [-xxyy_avr;
     -xxzz_avr;
     -xxx_avr;
     -xxy_avr;
     -xxz_avr;
     -xx_avr];

resoult = inv(A0)*b;
x00 = -resoult(3)/2;                  %��ϳ���x����
y00 = -resoult(4)/(2*resoult(1));     %��ϳ���y����
z00 = -resoult(5)/(2*resoult(2));     %��ϳ���z����

AA = sqrt(x00*x00 + resoult(1)*y00*y00 + resoult(2)*z00*z00 - resoult(6));   % ��ϳ���x�����ϵ���뾶
BB = AA/sqrt(abs(resoult(1)));                                                    % ��ϳ���y�����ϵ���뾶
CC = AA/sqrt(abs(resoult(2)));                                                    % ��ϳ���z�����ϵ���뾶
 for m=1:num_points
    MagDataOut(m,1)=(x(m,1)-x00)/AA;
    MagDataOut(m,2)=(y(m,1)-y00)/BB;
    MagDataOut(m,3)=(z(m,1)-z00)/CC;
 end

end