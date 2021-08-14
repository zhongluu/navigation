%%%%%%%%%%%%%������ģ��%%%%%%%%%%%%%%%
 % Declaration:
 %  Copyright(c) 2021-2025, Navigation Research Center All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%%   the mode of gyroscope
%% ���������Wibb--���������
%% ���������atti--��������������򣨵�λ���ȣ�
%            atti_rate--������ʡ��������ʡ��������ʣ���λ����/�룩
%            t--����ʱ�䣨��λ���룩
%            T--���沽������λ���룩
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Wibb ] = real_gyroscope(atti,atti_rate)
 %% %%%%%%%%%%%%%������йصĲ���%%%%%%%%%%%%%%%%%%%%%%
%   Re=6378137.0;        %����뾶     ����λ���ף� 
%   f=1/298.257;         %�������Բ��
%   Wie=7.292115147e-5;  %������ת���ٶȣ���λ������/�룩
%   g=9.7803698;         %�������ٶ�    ����λ����/��/�룩
  %% ������λ��
%   long=posi(1,1)*pi/180.0;lati=posi(2,1)*pi/180.0;heig=posi(3,1);
  %% �������ʰ뾶���
%   Rm=Re*(1-2*f+3*f*sin(lati)*sin(lati));
%   Rn=Re*(1+f*sin(lati)*sin(lati));
  %% ��̬�Ǻ���̬������
  roll=atti(1,1)*pi/180.0;pitch=atti(2,1)*pi/180.0;head=atti(3,1)*pi/180.0;
  droll=atti_rate(1,1)*pi/180.0;dpitch=atti_rate(2,1)*pi/180.0; dhead=atti_rate(3,1)*pi/180.0;
  %% ����ϵN-->B
%   Cbn=[cos(roll)*cos(head)+sin(roll)*sin(pitch)*sin(head), -cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head), -sin(roll)*cos(pitch);
%        cos(pitch)*sin(head),                               cos(pitch)*cos(head),                                sin(pitch);
%        sin(roll)*cos(head)-cos(roll)*sin(pitch)*sin(head), -sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head), cos(roll)*cos(pitch)];
  %% ŷ���Ǳ任����
  Eluer_M=[cos(roll), 0, sin(roll)*cos(pitch);
           0,         1, -sin(pitch);
           sin(roll), 0, -cos(pitch)*cos(roll)];   
 
  %% %%%%%%%%%���������%%%%%%%%%%%%
  %% ��ز�������
  Wnbb=Eluer_M*[dpitch;droll;dhead];
%   veloN=Cbn'*veloB;
%   Ve=veloN(1,1);Vn=veloN(2,1);Vu=veloN(3,1); %Ve:���ȷ��� Vn:γ�ȷ��� Vu:�߶ȷ���
%   Wenn=[-Vn/(Rm+heig); Ve/(Rn+heig);  Ve/(Rn+heig)*tan(lati)];%����ϵ��Ե���ϵ��ת���������ڵ���ϵ�ϵ�ͶӰ
%   Wien=[0;             Wie*cos(lati); Wie*sin(lati)];
  %% ��ֵ�ļ���
  Wibb=Wnbb; %��λ������/��  ������ʵ���������ʾ��װ�ڻ���ϵ������������ڹ��Կռ�Ľ������ڻ���ϵ�ϵ�ͶӰ
  Wibb=Wibb*180.0/pi;        %��λ����/��
end




