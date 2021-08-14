function [ Pq,sq,Neff,RegKG,WibbAWR ] = initalparticleswithoutKLD(t,T,sq,Wibb,Fb,mag,RegKG,WibbAWR,noise_wARW,noise_gyro,noise_accel,noise_mag)
 %  Description:
 %   quaternion particle filter algorithm
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
g=9.7803698;
refer_mag=[0.48e-3;0.0;0.48e-3].*1;
refer_mag=[-0.009e-3;0.089e-3;0.043e-3];%东北天
  if t<30*T
      Qn=30000;
  elseif t<60*T
      Qn=20000;
  elseif t<80*T
      Qn=15000;      
  elseif t<100*T
      Qn=10000;
  elseif t<120*T
      Qn=6000;
  elseif t<150*T
      Qn=3000;      
  elseif t<180*T
      Qn=1500;
  elseif t<200*T
      Qn=800;
  elseif t<220*T
      Qn=400;
  elseif t<240*T
      Qn=100;
  elseif t<260*T
      Qn=30;
  elseif t<280*T
      Qn=5;  
  else
      Qn=1;
  end
for m=1:length(sq)
    WibbAWR(:,m)=WibbAWR(:,m)+noise_wARW*T*randn(3,1);
    [sq(:,m)]=quaternion_updata(T,Wibb-WibbAWR(:,m)*180/pi-rad2deg(noise_gyro)*Qn*randn(3,1),sq(:,m));%
end

if t<5*T
    R=800000;
elseif t<10*T
    R=700000;
elseif t<15*T
    R=600000;
elseif t<20*T
    R=500000;
elseif t<25*T
    R=400000;
elseif t<30*T
    R=300000;
elseif t<35*T
    R=200000;    
elseif t<40*T
    R=100000;
elseif t<45*T
    R=90000;
elseif t<60*T
    R=80000;
elseif t<75*T
    R=70000;
elseif t<90*T
    R=50000;
elseif t<100*T
    R=30000;
elseif t<110*T
    R=20000;
elseif t<115*T
    R=10000;
elseif t<125*T
    R=8000;
elseif t<130*T
    R=6000;
elseif t<135*T
    R=4000;    
elseif t<140*T
    R=2000;
elseif t<150*T
    R=1500;
elseif t<155*T
    R=1200;    
elseif t<160*T
    R=1000;   
elseif t<175*T
    R=800;  
elseif t<180*T
    R=600;  
elseif t<185*T
    R=500;  
elseif t<190*T
    R=400; 
elseif t<195*T
    R=200; 
elseif t<200*T
    R=100; 
elseif t<205*T
    R=80; 
elseif t<210*T
    R=60; 
elseif t<225*T
    R=50; 
elseif t<235*T
    R=40; 
elseif t<240*T
    R=30; 
elseif t<245*T
    R=20; 
elseif t<250*T
    R=10; 
elseif t<260*T
    R=8; 
elseif t<265*T
    R=6; 
elseif t<275*T
    R=4; 
elseif t<280*T
    R=2; 
else
    R=1;
end
for m=1:length(sq)
     Cnb=[sq(2,m)^2+sq(1,m)^2-sq(4,m)^2-sq(3,m)^2, 2*(sq(2,m)*sq(3,m)+sq(1,m)*sq(4,m)), 2*(sq(2,m)*sq(4,m)-sq(1,m)*sq(3,m));
       2*(sq(2,m)*sq(3,m)-sq(1,m)*sq(4,m)), sq(3,m)^2-sq(4,m)^2+sq(1,m)^2-sq(2,m)^2,  2*(sq(3,m)*sq(4,m)+sq(1,m)*sq(2,m));
       2*(sq(2,m)*sq(4,m)+sq(1,m)*sq(3,m)), 2*(sq(3,m)*sq(4,m)-sq(1,m)*sq(2,m)), sq(4,m)^2-sq(3,m)^2-sq(2,m)^2+sq(1,m)^2]; 
      KG(m)=RegKG(m)*mvnpdf([Cnb,zeros(3);zeros(3),Cnb]*[refer_mag;[0 0 -g]'],[mag;Fb],R.*[noise_mag.^2*eye(3),zeros(3);zeros(3),noise_accel.^2*eye(3)]);
end

%% 归一化权值
KGall=sum(KG);
if KGall==0
    RegKG=(1/length(KG))*ones(1,length(KG));
    [sq]=MGC_Generation( Eluer_to_quaternion([0 0 0]'),2,length(sq)/2,5000 );
else
    RegKG=KG./KGall;
end

%% 状态更新
Pq = qwaboE(sq,RegKG);% 求估计四元数均值
Neff=1/sum(RegKG.^2);
[ sq,WibbAWR,RegKG ] = resampler(sq,WibbAWR,RegKG,2*length(sq)/3);%,
end

