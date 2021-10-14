 %  Description:
 %   Quaternion weighted average based on eigenvector
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
%% 输入参数：qs--qsample采样四元素---4行n列
%             wis--采样四元素的归一化权值--1行n列
%% 输出参数：q--估计的四元素
%% 中间参数：M=sigma（wi*qi*qi'）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ q ] = qwaboE(qs,wis)
%% 计算采样的总数
[r,c]=size(qs);
n=c;
%%
M=zeros(r,r);
for i=1:n
    qi=qs(:,i);
    M=M+wis(n)*(qi)*(qi');
end
[V,D]=eig(M);
maxofD=max(D);
maxvalue=max(maxofD);
for i=1:c
    if(maxofD(i)==maxvalue)
        index=i;
        break
    end
end
q=V(:,index);
end

