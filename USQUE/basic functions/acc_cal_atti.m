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
%% ���������accN--accelerometer_only_g�������--ֻ���������еļ��ٶȼ����ֵ--����������3��1�У�
%% ���������attiN0--���Ƶ���̬�ǣ�3��1�У�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ attiNO ] = acc_cal_atti( accN )
%%
rad_to_deg=180/pi;
R=sqrt( accN(1,1)^2 + accN(2,1)^2 + accN(3,1)^2);%�����������ٶ�
AX=accN(1,1)/R;AY=accN(2,1)/R;AZ=accN(3,1)/R;%X����Y��ͷ��Z����
attiNO(1,1) = atan(AX/sqrt(AY^2+AZ^2))*rad_to_deg;%���
attiNO(2,1) = -atan(AY/sqrt(AX^2+AZ^2))*rad_to_deg;%����
attiNO(3,1) = atan(AZ/sqrt(AX^2+AY^2))*rad_to_deg;%����
end

