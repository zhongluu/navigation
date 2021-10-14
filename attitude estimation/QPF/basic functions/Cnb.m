 %  Description:
 %   attitude tranform matrix
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
 
function  ABC=Cnb(q0,q1,q2,q3)
   ABC=[q0^2+q1^2-q2^2-q3^2,  2*q1*q2+2*q0*q3,     2*q1*q3-2*q0*q2;
        2*q1*q2-2*q0*q3,      q0^2-q1^2+q2^2-q3^2, 2*q2*q3+2*q0*q1;
        2*q1*q3+2*q0*q2,      2*q3*q2-2*q0*q1,     q0^2-q1^2-q2^2+q3^2];   %77£ºµ¼º½p95£¨2.102£©

end