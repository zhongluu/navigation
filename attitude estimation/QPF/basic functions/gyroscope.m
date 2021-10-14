%%%%%%%%%%%%%陀螺仪模型%%%%%%%%%%%%%%%
 % Declaration:
 %  Copyright(c) 2021-2025, Navigation Research Center All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%%   the mode of gyroscope
%% 输出参数：Wibb--输出角速率
%            Gyro_r--陀螺一阶马尔可夫过程（弧度/秒）
%% 输入参数：atti--横滚、俯仰、航向（单位：度）
%            atti_rate--横滚速率、俯仰速率、航向速率（单位：度/秒）
%            acceB--飞机运动加速度——X右翼、Y机头、Z天向（单位：米/秒/秒）
%            posi--航迹发生器位置经度、纬度、高度（单位：度、度、米）
%            t--仿真时间（单位：秒）
%            T--仿真步长（单位：秒）
%            Gyro_r--陀螺一阶马尔可夫过程（弧度/秒）
%% 中间参数：Gyro_b--陀螺随机常数（弧度/秒）
%            Gyro_wg--陀螺白噪声
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Wibb,Gyro_r ] = gyroscope(t,T,atti,atti_rate,Gyro_r,noise_wARW,noise_gyro)
%% 姿态角和姿态角速率
roll=atti(1,1)*pi/180.0;pitch=atti(2,1)*pi/180.0;head=atti(3,1)*pi/180.0;
droll=atti_rate(1,1)*pi/180.0;dpitch=atti_rate(2,1)*pi/180.0; dhead=atti_rate(3,1)*pi/180.0;
%% 欧拉角变换矩阵
Eluer_M=[cos(roll), 0, sin(roll)*cos(pitch);
    0,         1, -sin(pitch);
    sin(roll), 0, -cos(pitch)*cos(roll)];
%% %%%%%%%%%陀螺仪输出%%%%%%%%%%%%
%% 相关参数计算
Wnbb=Eluer_M*[dpitch;droll;dhead];
%% 真值的计算
Wibb=Wnbb; %单位：弧度/秒  陀螺仪实际输出，表示安装在机体系的陀螺仪相对于惯性空间的角速率在机体系上的投影
Wibb=Wibb*180.0/pi;        %单位：度/秒
%% %%%%%%%%%加入误差%%%%%%%%%%%%%%
%% 参数设置
deg_rad=0.01745329252e0;% Transfer from angle degree to rad pi/180

Tg=3600.0; %陀螺一阶马尔可夫过程相关时间
Gyro_b=zeros(3,1);    %陀螺随机常数
Gyro_wg=zeros(3,1);   %陀螺白噪声
%% 误差产生
if( t==0 )
    Gyro_r=Gyro_r+noise_wARW*T*randn(3,1);
    Gyro_wg= noise_gyro*randn(3,1);
else
    Gyro_r=Gyro_r+noise_wARW*T*randn(3,1);%Gyro_r为陀螺仪漂移
    Gyro_wg=noise_gyro*randn(3,1);%0.3162e-6*randn(3,1);
end
%% 真实含误差输出
Wibb=Wibb+Gyro_wg/deg_rad+Gyro_r/deg_rad+Gyro_b/deg_rad;
end




