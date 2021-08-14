%%%%%%%%%%%%%陀螺仪模型%%%%%%%%%%%%%%%
 % Declaration:
 %  Copyright(c) 2021-2025, Navigation Research Center All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%%   the mode of gyroscope
%% 输出参数：Wibb--输出角速率
%% 输入参数：atti--横滚、俯仰、航向（单位：度）
%            atti_rate--横滚速率、俯仰速率、航向速率（单位：度/秒）
%            t--仿真时间（单位：秒）
%            T--仿真步长（单位：秒）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Wibb ] = real_gyroscope(atti,atti_rate)
 %% %%%%%%%%%%%%%与地球有关的参数%%%%%%%%%%%%%%%%%%%%%%
%   Re=6378137.0;        %地球半径     （单位：米） 
%   f=1/298.257;         %地球的椭圆率
%   Wie=7.292115147e-5;  %地球自转角速度（单位：弧度/秒）
%   g=9.7803698;         %重力加速度    （单位：米/秒/秒）
  %% 飞行器位置
%   long=posi(1,1)*pi/180.0;lati=posi(2,1)*pi/180.0;heig=posi(3,1);
  %% 地球曲率半径求解
%   Rm=Re*(1-2*f+3*f*sin(lati)*sin(lati));
%   Rn=Re*(1+f*sin(lati)*sin(lati));
  %% 姿态角和姿态角速率
  roll=atti(1,1)*pi/180.0;pitch=atti(2,1)*pi/180.0;head=atti(3,1)*pi/180.0;
  droll=atti_rate(1,1)*pi/180.0;dpitch=atti_rate(2,1)*pi/180.0; dhead=atti_rate(3,1)*pi/180.0;
  %% 坐标系N-->B
%   Cbn=[cos(roll)*cos(head)+sin(roll)*sin(pitch)*sin(head), -cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head), -sin(roll)*cos(pitch);
%        cos(pitch)*sin(head),                               cos(pitch)*cos(head),                                sin(pitch);
%        sin(roll)*cos(head)-cos(roll)*sin(pitch)*sin(head), -sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head), cos(roll)*cos(pitch)];
  %% 欧拉角变换矩阵
  Eluer_M=[cos(roll), 0, sin(roll)*cos(pitch);
           0,         1, -sin(pitch);
           sin(roll), 0, -cos(pitch)*cos(roll)];   
 
  %% %%%%%%%%%陀螺仪输出%%%%%%%%%%%%
  %% 相关参数计算
  Wnbb=Eluer_M*[dpitch;droll;dhead];
%   veloN=Cbn'*veloB;
%   Ve=veloN(1,1);Vn=veloN(2,1);Vu=veloN(3,1); %Ve:经度方向 Vn:纬度方向 Vu:高度方向
%   Wenn=[-Vn/(Rm+heig); Ve/(Rn+heig);  Ve/(Rn+heig)*tan(lati)];%导航系相对地球系的转动角速率在导航系上的投影
%   Wien=[0;             Wie*cos(lati); Wie*sin(lati)];
  %% 真值的计算
  Wibb=Wnbb; %单位：弧度/秒  陀螺仪实际输出，表示安装在机体系的陀螺仪相对于惯性空间的角速率在机体系上的投影
  Wibb=Wibb*180.0/pi;        %单位：度/秒
end




