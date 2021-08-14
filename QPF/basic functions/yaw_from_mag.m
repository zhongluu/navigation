function [ yaw ] = yaw_from_mag( roll,pitch,mag )
 %  Description:
 %   yaw drived from magnet
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
roll=roll*pi/180;pitch=pitch*pi/180;
mx=mag(1)*cos(roll)+mag(3)*sin(roll);
my=mag(1)*sin(pitch)*sin(roll)+mag(2)*cos(pitch)-mag(3)*cos(roll)*sin(pitch);
yaw=-atan(mx/my)*180/pi;
end

