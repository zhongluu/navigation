 %  Description:
 %   Anti-Matrix
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

function [ A ] = qAntiMatrix( q )
[r,~]=size(q);
if(r==4)
A=[q(1) -q(2) -q(3) -q(4);
    q(2)  q(1) -q(4)  q(3);
    q(3)  q(4)  q(1) -q(2);
    q(4) -q(3)  q(2)  q(1)];
elseif(r==3)
    A=[0,-q(3),q(2);
       q(3),0,-q(1);
      -q(2),q(1),0];
end
end

