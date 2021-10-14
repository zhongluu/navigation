%%%%%%%%%%%%%������ģ��%%%%%%%%%%%%%%%
 % Declaration:
 %  Copyright(c) 2021-2025, Navigation Research Center All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%%   the mode of gyroscope
%% ���������Wibb--���������
%            Gyro_r--����һ������ɷ���̣�����/�룩
%% ���������atti--��������������򣨵�λ���ȣ�
%            atti_rate--������ʡ��������ʡ��������ʣ���λ����/�룩
%            acceB--�ɻ��˶����ٶȡ���X����Y��ͷ��Z���򣨵�λ����/��/�룩
%            posi--����������λ�þ��ȡ�γ�ȡ��߶ȣ���λ���ȡ��ȡ��ף�
%            t--����ʱ�䣨��λ���룩
%            T--���沽������λ���룩
%            Gyro_r--����һ������ɷ���̣�����/�룩
%% �м������Gyro_b--�����������������/�룩
%            Gyro_wg--���ݰ�����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Wibb,Gyro_r ] = gyroscope(t,T,atti,atti_rate,Gyro_r,noise_wARW,noise_gyro)
%% ��̬�Ǻ���̬������
roll=atti(1,1)*pi/180.0;pitch=atti(2,1)*pi/180.0;head=atti(3,1)*pi/180.0;
droll=atti_rate(1,1)*pi/180.0;dpitch=atti_rate(2,1)*pi/180.0; dhead=atti_rate(3,1)*pi/180.0;
%% ŷ���Ǳ任����
Eluer_M=[cos(roll), 0, sin(roll)*cos(pitch);
    0,         1, -sin(pitch);
    sin(roll), 0, -cos(pitch)*cos(roll)];
%% %%%%%%%%%���������%%%%%%%%%%%%
%% ��ز�������
Wnbb=Eluer_M*[dpitch;droll;dhead];
%% ��ֵ�ļ���
Wibb=Wnbb; %��λ������/��  ������ʵ���������ʾ��װ�ڻ���ϵ������������ڹ��Կռ�Ľ������ڻ���ϵ�ϵ�ͶӰ
Wibb=Wibb*180.0/pi;        %��λ����/��
%% %%%%%%%%%�������%%%%%%%%%%%%%%
%% ��������
deg_rad=0.01745329252e0;% Transfer from angle degree to rad pi/180

Tg=3600.0; %����һ������ɷ�������ʱ��
Gyro_b=zeros(3,1);    %�����������
Gyro_wg=zeros(3,1);   %���ݰ�����
%% ������
if( t==0 )
    Gyro_r=Gyro_r+noise_wARW*T*randn(3,1);
    Gyro_wg= noise_gyro*randn(3,1);
else
    Gyro_r=Gyro_r+noise_wARW*T*randn(3,1);%Gyro_rΪ������Ư��
    Gyro_wg=noise_gyro*randn(3,1);%0.3162e-6*randn(3,1);
end
%% ��ʵ��������
Wibb=Wibb+Gyro_wg/deg_rad+Gyro_r/deg_rad+Gyro_b/deg_rad;
end




