%% EKF_EULER
 %  Description:
 %   EKF_EULER
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
%zhong yuly
%2019-12-5
%reference:1刘建业--导航系统理论及应用
%          2Anton J.Haug--贝叶斯估计与跟踪实用指南
%          3朱志宇--粒子滤波及应用
%notation:tex
%%
%inital Patt can set as (0.5 deg)^2
%inital Pbias can set as (0.2deg/hr)^2
function [ X,Pxx ] = EKFEuler( T,X,omega,obs,NofGyro,NofBI,Nofacc,Nofmag,Pxx)
%Input:X-->state=[pitch,roll,yaw,beta]'--6-dimension
%      T-->time interval
%      omega-->the output of gyro (deg/s)
%      obs-->observed vector:combined by the output of acc. and mag.
%      Nofwibb-->the deviation of noise of gyro (deg/s)
%      NofBI-->the devaition of bias instability noise of gyro
%      Nofacc-->the devaition of noise of accelerometer
%      Nofmag-->the devaition of noise of magnetmeter
%      P-->state covariance( priori )
%% required parameter
dpitch=omega(1,1)*pi/180.0;droll=omega(2,1)*pi/180.0; dhead=omega(3,1)*pi/180.0;
theta=X(1);gamma=X(2);

G=(1/cos(theta))*[cos(theta)*cos(gamma)          0         sin(gamma)*cos(theta);...
                     sin(gamma)*sin(theta)   cos(theta)    -cos(gamma)*sin(theta);...
                         sin(gamma)              0             -cos(gamma)];
tao=[-T*G    ,zeros(3);...
     zeros(3),T*eye(3) ];%
W=[NofGyro*eye(3)  ,zeros(3);  ...
    zeros(3)       ,NofBI*eye(3)];
Q=(tao*W)*(tao*W)';
R=[Nofmag^2*eye(3),zeros(3);...
    zeros(3)      ,Nofacc^2*eye(3)];
%% start
Feuler=eye(3)+T*[0                                                                            -sin(gamma)*dpitch+cos(gamma)*dhead                      0;...
          (sin(gamma)*dpitch+cos(gamma)*dhead)*sec(theta)^2                                   tan(theta)*(cos(gamma)*dpitch+sin(gamma)*dhead)          0;...
    (sin(gamma)*dpitch*tan(theta)/cos(theta))-(cos(gamma)*dhead*tan(theta)/cos(theta))       cos(gamma)/cos(theta)*dpitch+sin(gamma)/cos(theta)*dhead  0 ];%
% FgyroBI=(-T)*[cos(gamma)            0     sin(gamma);...
%               sin(gamma)*            1    -cos(gamma);...
%               sin(gamma)/cos(theta) 0    -cos(gamma)/cos(theta)];
F=[  Feuler ,(-T)*G;...
    zeros(3),eye(3) ];
%% prediction
X(1:3)=X(1:3)+(T*G)*[dpitch,droll,dhead]';
% X=F*X+tao*[NofGyro*randn(3,1);NofBI*randn(3,1)];
Pxx=F*Pxx*F'+Q;
%% correction
ref_g=9.7803698;
ref_mag=0.48e-3;
ref=[[ref_mag,0.0,ref_mag],[0 0 ref_g]]';
theta=X(1);gamma=X(2);psi=X(3);

Cnb=[cos(gamma)*cos(psi)+sin(gamma)*sin(theta)*sin(psi), -cos(gamma)*sin(psi)+sin(gamma)*sin(theta)*cos(psi), -sin(gamma)*cos(theta);...
     cos(theta)*sin(psi),                               cos(theta)*cos(psi),                                sin(theta);...
     sin(gamma)*cos(psi)-cos(gamma)*sin(theta)*sin(psi), -sin(gamma)*sin(psi)-cos(gamma)*sin(theta)*cos(psi), cos(gamma)*cos(theta)];


% Cnb=[cos(gamma)*cos(psi)+sin(gamma)*sin(theta)*sin(psi) -cos(gamma)*sin(psi)+sin(gamma)*sin(theta)*cos(psi) -sin(gamma)*cos(theta);...
%                 cos(theta)*sin(psi)                                  cos(theta)*cos(psi)                             sin(theta);...
%      sin(gamma)*cos(psi)-cos(gamma)*sin(theta)*sin(psi) -sin(gamma)*sin(psi)-cos(gamma)*sin(theta)*cos(psi)  cos(gamma)*cos(theta)];
Y=[Cnb      ,zeros(3);...
    zeros(3),Cnb    ]*ref;
H=[    ref_mag*[sin(gamma)*cos(theta)*sin(psi)+sin(gamma)*sin(theta) -sin(gamma)*cos(psi)+cos(gamma)*sin(theta)*sin(psi)-cos(gamma)*cos(theta)  cos(gamma)*sin(psi)+sin(gamma)*sin(theta)*cos(psi);...
               -sin(theta)*sin(psi)+cos(theta)                                        0                                                         cos(theta)*cos(psi);...
               -cos(gamma)*cos(theta)*sin(psi)-cos(gamma)*sin(theta) cos(gamma)*cos(psi)+sin(gamma)*sin(theta)*sin(psi)-sin(gamma)*cos(theta)  -sin(gamma)*sin(psi)-cos(gamma)*cos(psi)*sin(theta)]
    ref_g*[              sin(gamma)*sin(theta)                         -cos(gamma)*cos(theta)                                                                        0;...
                             cos(theta)                                        0                                                                                 0;...
                         -cos(gamma)*sin(theta)                                 -sin(gamma)*cos(theta)                                                               0];];
H=[H,zeros(6,3)];
Pyy=H*Pxx*H'+R;
Pxy=Pxx*H';
K=Pxy/Pyy;
X=X+K*(obs-Y);
Pxx=Pxx-K*Pyy*K';
end

