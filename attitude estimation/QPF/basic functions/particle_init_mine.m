 %  Description:
 %   quaternion particle initial
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

function [q]=particle_init_mine(b0,r0,M)
q=zeros(4,M);
noise_b0=2.8935e-4;
for i=1:M
    b=b0+noise_b0*randn(3,1);
    theta=0.5*acos(b'*r0/norm(b)*norm(r0));
    e=cross(b,r0)./norm(cross(b,r0));
    beta=2*pi*rand;
    q(:,i)=qAntiMatrix([cos(beta);sin(beta).*b0./norm(b0)])*[cos(theta);sin(theta)*e];
end
end