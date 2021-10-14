 %  Description:
 %    Jacobian transform
 %  Reference:
 %  Crassidis, John L.Markley, F. Landis. Unscented Filtering for Spacecraft Attitude Estimation
 %  Oshman, Yaakov Carmi, Avishy.  Attitude Estimation from Vector Observations Using a Genetic-Algorithm-Embedded
 %                                     Quaternion Particle Filter
 %  Zhaihe ZHOU, Yulu ZHONG, Chuanwei ZENG, and Xiangrui TIAN. Attitude Estimation Using Parallel Quaternion Particle 
 %                                                                      Filter Based on New Quaternion Distribution
 % Declaration:
 %  Copyright(c) 2021-2025, by Yulu Zhong, Chuanwei Zeng, All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%%
%% 输入参数：uq--样本均值
%           sigmaq--样本方差
%% 输出参数：G--雅克比矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ G ] = Jacobian( uq )
G=zeros(3,4);
G(3,1)=(-(uq(3)+uq(2))/((uq(3)+uq(2))^2+(uq(4)+uq(1))^2))+((uq(3)-uq(2))/((uq(3)-uq(2))^2+(uq(4)-uq(1))^2));%ROLL
G(3,2)=( (uq(4)+uq(1))/((uq(3)+uq(2))^2+(uq(4)+uq(1))^2))-((uq(4)-uq(1))/((uq(3)-uq(2))^2+(uq(4)-uq(1))^2));
G(3,3)=( (uq(4)+uq(1))/((uq(3)+uq(2))^2+(uq(4)+uq(1))^2))+((uq(4)-uq(1))/((uq(3)-uq(2))^2+(uq(4)-uq(1))^2));
G(3,4)=(-(uq(3)+uq(2))/((uq(3)+uq(2))^2+(uq(4)+uq(1))^2))-((uq(3)-uq(2))/((uq(3)-uq(2))^2+(uq(4)-uq(1))^2));
G(2,1)=(2*uq(4)/sqrt(1-4*(uq(2)*uq(3)+uq(1)*uq(4))^2));%PITCH
G(2,2)=(2*uq(3)/sqrt(1-4*(uq(2)*uq(3)+uq(1)*uq(4))^2));
G(2,3)=(2*uq(2)/sqrt(1-4*(uq(2)*uq(3)+uq(1)*uq(4))^2));
G(2,4)=(2*uq(1)/sqrt(1-4*(uq(2)*uq(3)+uq(1)*uq(4))^2));
G(1,1)=(-(uq(3)+uq(2))/((uq(3)+uq(2))^2+(uq(4)+uq(1))^2))-((uq(3)-uq(2))/((uq(3)-uq(2))^2+(uq(4)-uq(1))^2));%YAW
G(1,2)=( (uq(4)+uq(1))/((uq(3)+uq(2))^2+(uq(4)+uq(1))^2))+((uq(4)-uq(1))/((uq(3)-uq(2))^2+(uq(4)-uq(1))^2));
G(1,3)=( (uq(4)+uq(1))/((uq(3)+uq(2))^2+(uq(4)+uq(1))^2))-((uq(4)-uq(1))/((uq(3)-uq(2))^2+(uq(4)-uq(1))^2));
G(1,4)=(-(uq(3)+uq(2))/((uq(3)+uq(2))^2+(uq(4)+uq(1))^2))+((uq(3)-uq(2))/((uq(3)-uq(2))^2+(uq(4)-uq(1))^2));
end

