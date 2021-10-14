%% USQUE
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
%zhong yulu
%reference:Unscented Filtering For Spacecraft Attitude Estimation-John
%L.Crassidis-AIAA
%2019-11-25
%notation similar tex
%% prior information
%inital state must set as [zeros(1,3) beta0]'
%inital Patt can set as (0.5 deg)^2
%inital Pbias can set as (0.2deg/hr)^2
%
function [Opt_q,X,P] = USQUE(T,X,q,omega,obs,Nofwibb,NofBI,Nofacc,Nofmag,P)
%Input:X-->state=[delta_p,beta]' has 6-dimension 
%      q-->previous quaternion
%      t-->time stamp
%      T-->time interval
%      omega-->the output of gyro (deg/s)
%      obs-->observed vector:combined by the output of acc. and mag.
%      Nofwibb-->the deviation of noise of gyro (deg/s)
%      NofBI-->the devaition of bias instability noise of gyro
%      Nofacc-->the devaition of noise of accelerometer
%      Nofmag-->the devaition of noise of magnetmeter
%      P-->state covariance( priori )
%Output:X--> the updated state
%       Opt_q-->the optimal quaternion 
%       P-->state covariance( posteriori )
[~,col]=size(X);
if col>1
    error('state must be column vector');
end
%% required parameter
a=1;
f=2*(a+1);
Q=(T/2)*[(Nofwibb^2-(1/6)*NofBI^2*T^2)*eye(3), zeros(3);...
         zeros(3)            ,          NofBI^2*eye(3)];
R=[Nofmag^2*eye(3),zeros(3);...
    zeros(3)      ,Nofacc^2*eye(3)];
n=6;
lambda=1;
%% start
sigma=[zeros(size(X)),-1*chol((n+lambda)*(P+Q)),chol((n+lambda)*(P+Q))];
sigmadots=length(sigma);
aleph=repmat(X,1,sigmadots)+sigma;
%convert MRPs to erro quaternion
delta_q=zeros(4,sigmadots);
for i=1:sigmadots
   scalar_delta_q=(-a*norm(aleph(1:3,i))^2+f*sqrt(f^2+(1-a^2)*norm(aleph(1:3,i))^2))/(f^2+norm(aleph(1:3,i))^2);
   vector_delta_q=(1/f)*(a+scalar_delta_q).*aleph(1:3,i);
   delta_q(:,i)=[scalar_delta_q;vector_delta_q];
%    delta_q(:,i)=delta_q(:,i)./norm(delta_q(:,i));
end
sq=zeros(4,sigmadots);
for i=1:sigmadots
    sq(:,i) =qAntiMatrix(delta_q(:,i))*q;
end
est_omega=repmat(omega,1,sigmadots)-aleph(4:6,:).*180/pi;
%% prediction
pred_q=zeros(4,sigmadots);
for i=1:sigmadots
    pred_q(:,i)=quaternion_updata(T,est_omega(:,i),sq(:,i));
end
pred_delta_q=zeros(4,sigmadots);
for i=1:sigmadots
    pred_delta_q(:,i)=qAntiMatrix(pred_q(:,i))*[pred_q(1,1);-pred_q(2:4,1)];
end
%convert error quaternion to MRPs
for i=1:sigmadots
    aleph(1:3)=f*(pred_delta_q(2:4)./(a+pred_delta_q(1)));
end
pred_mean=(1/(n+lambda))*(lambda*aleph(:,1)+(1/2)*sum(aleph(:,2:end),2));
pred_cov=(1/(n+lambda))*(lambda*(aleph(:,1)-pred_mean)*(aleph(:,1)-pred_mean)'...
    +(1/2)*((aleph(:,2:end)-repmat(pred_mean,1,sigmadots-1))*(aleph(:,2:end)-repmat(pred_mean,1,sigmadots-1))'))+Q;
%% correction
g=9.7803698; %gravity acceleration
ref=[[0.48e-3,0.0,0.48e-3],[0 0 -g]]';%combined by mag. and acc.
gamma=zeros(length(ref),sigmadots);
for i=1:sigmadots
    Cnb=[pred_q(2,i)^2+pred_q(1,i)^2-pred_q(4,i)^2-pred_q(3,i)^2, 2*(pred_q(2,i)*pred_q(3,i)+pred_q(1,i)*pred_q(4,i))     , 2*(pred_q(2,i)*pred_q(4,i)-pred_q(1,i)*pred_q(3,i));...
         2*(pred_q(2,i)*pred_q(3,i)-pred_q(1,i)*pred_q(4,i))    , pred_q(3,i)^2-pred_q(4,i)^2+pred_q(1,i)^2-pred_q(2,i)^2 , 2*(pred_q(3,i)*pred_q(4,i)+pred_q(1,i)*pred_q(2,i));...
         2*(pred_q(2,i)*pred_q(4,i)+pred_q(1,i)*pred_q(3,i))    , 2*(pred_q(3,i)*pred_q(4,i)-pred_q(1,i)*pred_q(2,i))     , pred_q(4,i)^2-pred_q(3,i)^2-pred_q(2,i)^2+pred_q(1,i)^2];
    gamma(:,i)=[Cnb,zeros(3);zeros(3),Cnb]*ref;
end
obs_mean=(1/(n+lambda))*(lambda*gamma(:,1)+(1/2)*sum(gamma(:,2:end),2));
obs_cov=(1/(n+lambda))*(lambda*(gamma(:,1)-obs_mean)*(gamma(:,1)-obs_mean)'...
    +(1/2)*((gamma(:,2:end)-repmat(obs_mean,1,sigmadots-1))*(gamma(:,2:end)-repmat(obs_mean,1,sigmadots-1))'));
innov_cov=obs_cov+R;
cross_corr_cov=(1/(n+lambda))*(lambda*(aleph(:,1)-pred_mean)*(gamma(:,1)-obs_mean)'...
    +(1/2)*((aleph(:,2:end)-repmat(pred_mean,1,sigmadots-1))*(gamma(:,2:end)-repmat(obs_mean,1,sigmadots-1))'));
KGain=cross_corr_cov/innov_cov;
X=X+KGain*(obs-obs_mean);
P=pred_cov-KGain*innov_cov*KGain';
%estimate the optimal quaternion
scalar_opt_delta_q=(-a*norm(X(1:3))^2+f*sqrt(f^2+(1-a^2)*norm(X(1:3))^2))/(f^2+norm(X(1:3))^2);
vector_delta_q=(1/f)*(a+scalar_opt_delta_q).*X(1:3);
opt_delta_q=[scalar_opt_delta_q ;vector_delta_q];
Opt_q=qAntiMatrix(opt_delta_q)*pred_q(:,1);
% Opt_q=Opt_q./norm(Opt_q);
%reset the MRPs
X(1:3)=zeros(3,1);
end

