function [ sigma_q ] = sigma_updata(T,Wibb,uq,sigma_q,nosiewibbvar)
 %  Description:
 %   sigma update
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
Wnbb=Wibb*pi/180.0;  %单位：弧度/秒
%%%%%%%%%
WnbbA=Wnbb*T;
WnbbA0=sqrt(WnbbA(1,1)^2+WnbbA(2,1)^2+WnbbA(3,1)^2);
WnbbX=[0,          -WnbbA(1,1), -WnbbA(2,1), -WnbbA(3,1);
       WnbbA(1,1),  0,           WnbbA(3,1), -WnbbA(2,1);
       WnbbA(2,1), -WnbbA(3,1),   0,          WnbbA(1,1);
       WnbbA(3,1),  WnbbA(2,1),  -WnbbA(1,1),   0         ];
c_q=cos(WnbbA0/2);
if( WnbbA0<=1.0e-15 ) d_q=0.5; else  d_q=sin(WnbbA0/2)/WnbbA0;  end
 A=c_q*eye(4)+d_q*WnbbX ;
sigma_q=A*sigma_q*A';
% Tg=qAntiMatrix(uq(2:4));
% Tg=[ -(uq(2:4))'
%     Tg+eye(3)*uq(1)];
% Tg=0.5*Tg*T;%计算系统噪声矩阵

% % sigma_q=( c_q*eye(4)+d_q*WnbbX )*sigma_q*( c_q*eye(4)+d_q*WnbbX )'+Tg*diag(nosiewibbvar)*Tg';%一般迭代算法 

end