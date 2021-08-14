 %  Description:
 %   Standardized Particle Filter
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
%% 输入参数：M--粒子数量
%            sq--采样粒子集
%            RegKG--粒子权值
%            Wibb--陀螺仪输出值
%            Fb--加速度计输出值
%            nosiewibbvar--陀螺仪方差
%% 输出参数：Pq--估计四元数
%            sq--重采样后的粒子集
%            RegKG--采样后权值
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Pq,sq,RegKG,WibbAWR,flag]= SPF_of_atti2( T,M,sq,RegKG,Wibb_last,Wibb,Fb,mag,WibbAWR,noise_wARW,noise_gyro,noise_mag, noise_accel)
g=9.7803698; 

refer_mag=[-0.009e-3;0.089e-3;0.043e-3];%东北天
correlation=[];correlation1=[];correlation2=[];
flag=0;
%% 时间更新
for m=1:M
    WibbAWR(:,m)=WibbAWR(:,m)+noise_wARW*T*randn(3,1);
    [sq(:,m)]=quaternion_updata(T,Wibb-rad2deg(WibbAWR(:,m))-rad2deg(noise_gyro)*randn(3,1),sq(:,m));%
end
%% 权值更新
for m=1:M

     Cnb=[sq(2,m)^2+sq(1,m)^2-sq(4,m)^2-sq(3,m)^2, 2*(sq(2,m)*sq(3,m)+sq(1,m)*sq(4,m)), 2*(sq(2,m)*sq(4,m)-sq(1,m)*sq(3,m));
       2*(sq(2,m)*sq(3,m)-sq(1,m)*sq(4,m)), sq(3,m)^2-sq(4,m)^2+sq(1,m)^2-sq(2,m)^2,  2*(sq(3,m)*sq(4,m)+sq(1,m)*sq(2,m));
       2*(sq(2,m)*sq(4,m)+sq(1,m)*sq(3,m)), 2*(sq(3,m)*sq(4,m)-sq(1,m)*sq(2,m)), sq(4,m)^2-sq(3,m)^2-sq(2,m)^2+sq(1,m)^2]; 

     Liangce=[Cnb,zeros(3);zeros(3),Cnb]*[refer_mag;[0 0 -g]'];
     Liangce2=[Cnb]*[refer_mag];
     Measurement=[mag;Fb];Measurement2=[mag];

     l(m)=mvnpdf(Liangce,Measurement,[noise_mag.^2*eye(3),zeros(3);zeros(3),noise_accel.^2*eye(3)]);
     KG(m)=RegKG(m)*l(m);

end
%% 归一化权值
KGall=sum(KG);
RegKG=KG./KGall;
%% 状态更新
Pq = qwaboE(sq,RegKG);% 求估计四元数均值
varPq = varofPq(sq,Pq,RegKG);% 求估计四元数方差
%% 重采样
[ sq,WibbAWR,RegKG ] = resampler(sq,WibbAWR,RegKG,2*M/3);

end

