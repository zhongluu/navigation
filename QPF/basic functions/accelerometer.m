%%%%%%%%%%%%���ٶȼ�ģ��%%%%%%%%%%%%
 % Declaration:
 %  Copyright(c) 2021-2025, Navigation Research Center All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%%    the mode of accelerometer
%% ���������Fb--����ϵ���ٶȼ���� ����λ����/��/�룩
%            Acc_r--���ٶ�һ������ɷ���̣�����/�룩
%% ���������T--���沽��
%            atti--��������������򣨵�λ���ȣ�
%            atti_rate--������ʡ��������ʡ��������ʣ���λ����/�룩
%            veloB--�ɻ��˶��ٶȡ���X����Y��ͷ��Z���򣨵�λ����/�룩
%            acceB--�ɻ��˶����ٶȡ���X����Y��ͷ��Z���򣨵�λ����/��/�룩
%            posi--����������λ�þ��ȡ�γ�ȡ��߶ȣ���λ���ȡ��ȡ��ף�
%            Acc_r--���ٶ�һ������ɷ���̣�����/�룩
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fb,Acc_r] = accelerometer(t,T,atti,noise_accel,Acc_r)
  %% %%%%%%%%%%%%%������йصĲ���%%%%%%%%%%%%%%%%%%%%%%
  g=9.7803698;         %�������ٶ�    ����λ����/��/�룩
  %% ��̬�Ǻ���̬������
  roll=atti(1,1)*pi/180.0;pitch=atti(2,1)*pi/180.0;head=atti(3,1)*pi/180.0;
  %% ����ϵN-->B
  Cbn=[cos(roll)*cos(head)+sin(roll)*sin(pitch)*sin(head), -cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head), -sin(roll)*cos(pitch);
       cos(pitch)*sin(head),                               cos(pitch)*cos(head),                                sin(pitch);
       sin(roll)*cos(head)-cos(roll)*sin(pitch)*sin(head), -sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head), cos(roll)*cos(pitch)];
  %% %%%%%%%%%%%%���ٶȼ����%%%%%%%%%%%%%%
  %% ��ֵ�ļ���
  Fn=-[0.0;0.0;1.0*g];  %��������
  Fb=Cbn*Fn;            %��λ����/��/��
  %% ��������
  deg_rad=0.01745329252e0;% Transfer from angle degree to rad 
  Da_bias=[1; 1; 1]*60e-6*g;  %���ٶ���ƫ��Ӧ������������еļ��ٶ���ƫ����һ�£�
  Ta=1800.0; %���ٶ�һ������ɷ�������ʱ��
  Acc_w=noise_accel*randn(3,1);
  %% �������
%   if( t==0 )
%       Acc_r=Da_bias.*randn(3,1); %���ٶ�һ������ɷ����1.0e-4g
%   else
%       Acc_wa=sqrt(2*T/Ta)*Da_bias.*randn(3,1);%���ٶ�һ������ɷ���̰�����
%       Acc_r=exp(-1.0*T/Ta)*Acc_r+Acc_wa; %���ٶ�һ������ɷ����
%   end
  %% ��ʵ��������
  Fb=Fb+Acc_w;%+Acc_r
end
