 %  Description:
 %   attitude eastime based on acclerometer
 %  Reference:
 %  Crassidis, John L.Markley, F. Landis. Unscented Filtering for Spacecraft Attitude Estimation
 %  Oshman, Yaakov Carmi, Avishy.  Attitude Estimation from Vector Observations Using a Genetic-Algorithm-Embedded
 %                                     Quaternion Particle Filter
 %  Zhaihe ZHOU, Yulu ZHONG, Chuanwei ZENG, and Xiangrui TIAN. Attitude Estimation Using Parallel Quaternion Particle 
 %                                                                      Filter Based on New Quaternion Distribution
 %  Yulu ZHONG. Research of Attitude Estimation Algorithm Based on Quaternion Fast Particle Filter
 % Declaration:
 %  Copyright(c) 2021-2025, by Yulu Zhong, Chuanwei Zeng, All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%% 
%% 输入参数：accN--accelerometer_only_g输出参数--只对重力敏感的加速度计输出值--理想条件（3行1列）
%% 输出参数：attiN0--估计的姿态角（3行1列）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ attiNO ] = acc_cal_atti( accN )
%%
rad_to_deg=180/pi;
R=sqrt( accN(1,1)^2 + accN(2,1)^2 + accN(3,1)^2);%计算重力加速度
AX=accN(1,1)/R;AY=accN(2,1)/R;AZ=accN(3,1)/R;%X右翼、Y机头、Z天向
attiNO(1,1) = atan(AX/sqrt(AY^2+AZ^2))*rad_to_deg;%横滚
attiNO(2,1) = -atan(AY/sqrt(AX^2+AZ^2))*rad_to_deg;%俯仰
attiNO(3,1) = atan(AZ/sqrt(AX^2+AY^2))*rad_to_deg;%航向
end

