%%%%%%%%%%%%加速度计模型%%%%%%%%%%%%
 % Declaration:
 %  Copyright(c) 2021-2025, Navigation Research Center All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%%    the mode of accelerometer
%% 输出参数：Fb--机体系加速度计输出 （单位：米/秒/秒）
%            Acc_r--加速度一阶马尔可夫过程（弧度/秒）
%% 输入参数：T--仿真步长
%            atti--横滚、俯仰、航向（单位：度）
%            atti_rate--横滚速率、俯仰速率、航向速率（单位：度/秒）
%            veloB--飞机运动速度——X右翼、Y机头、Z天向（单位：米/秒）
%            acceB--飞机运动加速度——X右翼、Y机头、Z天向（单位：米/秒/秒）
%            posi--航迹发生器位置经度、纬度、高度（单位：度、度、米）
%            Acc_r--加速度一阶马尔可夫过程（弧度/秒）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fb,Acc_r] = accelerometer(t,T,atti,noise_accel,Acc_r)
  %% %%%%%%%%%%%%%与地球有关的参数%%%%%%%%%%%%%%%%%%%%%%
  g=9.7803698;         %重力加速度    （单位：米/秒/秒）
  %% 姿态角和姿态角速率
  roll=atti(1,1)*pi/180.0;pitch=atti(2,1)*pi/180.0;head=atti(3,1)*pi/180.0;
  %% 坐标系N-->B
  Cbn=[cos(roll)*cos(head)+sin(roll)*sin(pitch)*sin(head), -cos(roll)*sin(head)+sin(roll)*sin(pitch)*cos(head), -sin(roll)*cos(pitch);
       cos(pitch)*sin(head),                               cos(pitch)*cos(head),                                sin(pitch);
       sin(roll)*cos(head)-cos(roll)*sin(pitch)*sin(head), -sin(roll)*sin(head)-cos(roll)*sin(pitch)*cos(head), cos(roll)*cos(pitch)];
  %% %%%%%%%%%%%%加速度计输出%%%%%%%%%%%%%%
  %% 真值的计算
  Fn=-[0.0;0.0;1.0*g];  %比力方程
  Fb=Cbn*Fn;            %单位：米/秒/秒
  %% 参数设置
  deg_rad=0.01745329252e0;% Transfer from angle degree to rad 
  Da_bias=[1; 1; 1]*60e-6*g;  %加速度零偏（应与非随机性误差中的加速度零偏保持一致）
  Ta=1800.0; %加速度一阶马尔可夫过程相关时间
  Acc_w=noise_accel*randn(3,1);
  %% 产生误差
%   if( t==0 )
%       Acc_r=Da_bias.*randn(3,1); %加速度一阶马尔可夫过程1.0e-4g
%   else
%       Acc_wa=sqrt(2*T/Ta)*Da_bias.*randn(3,1);%加速度一阶马尔可夫过程白噪声
%       Acc_r=exp(-1.0*T/Ta)*Acc_r+Acc_wa; %加速度一阶马尔可夫过程
%   end
  %% 真实含误差输出
  Fb=Fb+Acc_w;%+Acc_r
end
