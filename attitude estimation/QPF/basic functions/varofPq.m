 %  Description:
 %   compute the var of quaternion
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
%% 输入参数：qs--采样粒子--4行n列
%             Pq--估计四元数
%             wis--粒子权值--1行n列
%% 输出参数：varPq--估计协方差矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ varPq ] = varofPq(qs,Pq,wis)
[~,c]=size(qs);
n=c;
qs=[1,0,0,0;
    0,-1,0,0;
    0,0,-1,0;
    0,0,0,-1] * qs;
varPq=0;
for i=1:n
    varPq=varPq+wis(i)*(qAntiMatrix(Pq)*qs(:,i))*(qAntiMatrix(Pq)*qs(:,i))';
end
% for i=1:n
%     varPq=varPq+wis(i)*(qAntiMatrix(qs(:,i))*Pq)*(qAntiMatrix(qs(:,i))*Pq)';
% end
end

