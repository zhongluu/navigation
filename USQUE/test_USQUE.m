%%%%%%%%%%%%%USQUE%%%%%%%%%%%%%%%%%%%%%
 %  Description:
 %   quaternion particle filter algorithm
 %  Reference:
 %  Crassidis, John L.Markley, F. Landis. Unscented Filtering for Spacecraft Attitude Estimation
 %  Oshman, Yaakov Carmi, Avishy.  Attitude Estimation from Vector Observations Using a Genetic-Algorithm-Embedded
 %                                     Quaternion Particle Filter
 %  Zhaihe ZHOU, Yulu ZHONG, Chuanwei ZENG, and Xiangrui TIAN. Attitude Estimation Using Parallel Quaternion Particle 
 %                                                                      Filter Based on New Quaternion Distribution
 %  Yulu ZHONG. Research of Attitude Estimation Algorithm Based on Quaternion Fast Particle Filter
 % Declaration:
 %  Copyright(c) 2021-2025, by Yulu Zhong, Chuanwei Zeng, All rights reserved. 
 %    Nanjing University of Aeronautics and Astronautics, NanJing, P.R.China
 %  01/31/2020, 07/31/2025
%% test USQUE algorithm
%����¶
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;
addpath('basic functions');
% digits(1024);
%% ����ʱ������
n=0;
t=0;
T=0.01;
t_stop=50;
%% ����������
atti=zeros(3,1);     %��������������򣨵�λ���ȣ�
atti_rate=zeros(3,1);%������ʡ��������ʡ��������ʣ���λ����/�룩

atti(1,1)=0.0;
atti(2,1)=0.0;
atti(3,1)=30.0;      %��ʼ����Ƕȣ���λ���ȣ�
%% Ԥ����IMU���
Wibb=zeros(3,1);    %����ϵ���������   ����λ����/�룩
Fb=zeros(3,1);      %����ϵ���ٶȼ���� ����λ����/��/�룩
Gyro_r=zeros(3,1);  % ����һ������ɷ���̣�����/�룩
Acc_r =zeros(3,1);  % ���ٶ�һ������ɷ���̣���/��/�룩
real_Wibb=zeros(3,1);
%% sensors specification
noise_wARW=deg2rad(0.0035);% rad/s
noise_gyro=deg2rad(0.1);% rad/s
noise_accel=0.02;% m/s^2
noise_mag=(10e-7);% Gauss/s

%% sensors specification ��OpenIMU300BI��
noise_wARW=deg2rad(33e-4);% rad/s
noise_gyro=deg2rad(4.1667e-4);% rad/s
noise_accel=20e-6;% m/s^2
noise_mag=(10e-9);% Gauss/s
noise_wARW=deg2rad(35e-4);% rad/s
noise_gyro=deg2rad(0.08);% rad/s
noise_accel=0.004;% m/s^2
noise_mag=(25e-9);% Gauss/s
% %% sensors specification ��OpenIMU300ZI��
% noise_wARW=deg2rad(50e-4);
% noise_gyro=deg2rad(17e-4);% rad/s
% noise_accel=10e-6;% g
% noise_mag=25e-9;% T
%% sensors specification (MTI-G-710)
% noise_wARW=deg2rad(0.0095);
% noise_gyro=deg2rad(10/3600);% rad/s
% noise_accel=15e-6;% g/sqrt(Hz)
% noise_mag=50e-9;% T
%%  ��ʼ����Ԫ��
q= Eluer_to_quaternion( atti );
%% ��ʼ��UKF
ukf_q=Eluer_to_quaternion(atti+0.5*[-5 5 15]');%+2*randn(3,1)
% ukf_X=[zeros(1,3),rad2deg(noise_wARW)*randn(1,3)]';
ukf_X=[zeros(1,3),zeros(1,3)]';
ukf_p=[noise_gyro^2*eye(3),zeros(3);...
       zeros(3)    ,noise_wARW^2*eye(3)];
%% �������ݴ洢�ռ�
TraceData=[];
wibbdata=[];
real_qData=[];
ErroData=[];
%% ��ʼ����������̬��ز���
q_fromGyro=q;
atti_fromGyro=atti;
%% ��ʼ����������̬��ز���
real_q=q;
real_atti=atti;
%% ʵ������
acc_fb=[];gyro_Wibb=[];
data=xlsread('static measurement.xlsx');
acc_data=data(:,1:3)./4096;
% acc_data(:,1)=acc_data(:,1)-mean(acc_data(:,1));acc_data(:,2)=acc_data(:,2)-mean(acc_data(:,2));acc_data(:,3)=acc_data(:,3)-mean(acc_data(:,3));
acc_data=acc_data.*-9.7804;
gyro_data=data(:,4:6)*pi/180./16.4;%��λ������/�룩
% gyro_data(:,1)=gyro_data(:,1)-mean(gyro_data(:,1));gyro_data(:,2)=gyro_data(:,2)-mean(gyro_data(:,2));gyro_data(:,3)=gyro_data(:,3)-mean(gyro_data(:,3));
mag_data=data(1:12500,7:9).*0.3;
[mag_data_fit]=MagDataFitting(mag_data);
[PX,MODEL] = gmm(mag_data_fit(:,1), 1);
%% required parameter
a=1;
f=2*(a+1);
Q=(T/2)*[(noise_gyro^2-(1/6)*noise_wARW^2*T^2)*eye(3), zeros(3);...
         zeros(3)            ,          noise_wARW^2*eye(3)];
R=[noise_mag^2*eye(3),zeros(3);...
    zeros(3)      ,noise_accel^2*eye(3)];
n=6;
lambda=1;USQUE_t=zeros(1,100*t_stop+1);
%% ���濪ʼ
while t<=t_stop
    %% �������ݲ���
    [t,atti,atti_rate]=traceset2(t,T,atti,atti_rate);%�����������������й켣����
%     [t,atti,atti_rate]=traceset_test(t,T,atti,atti_rate);%�����������������й켣����
    [ Wibb,Gyro_r ] = gyroscope(t,T,atti,atti_rate,Gyro_r,noise_wARW,noise_gyro);%�������������������
    [ real_Wibb] = real_gyroscope(atti,atti_rate);%�������������������
    [real_q]=quaternion_updata(T,real_Wibb,real_q);%��������̬
    [Fb] =accelerometer(t,T,atti,noise_accel);%���������ٶȼ��������
    [mag]=magnet(atti,noise_mag);
    acc_fb=[acc_fb,Fb];
    gyro_Wibb=[gyro_Wibb,Wibb];
    %% ������ʵ��̬
    [Ture_q]=Eluer_to_quaternion( atti );
    
    %% ��������̬���
    [q_fromGyro]=quaternion_updata(T,Wibb,q_fromGyro);
    atti_fromGyro=atti_compute(q_fromGyro);
    [ atti_fromMeasure ] = acc_cal_atti( Fb );
    tmp_attiNO=deg2rad(atti_fromMeasure);
     atti_fromMeasure(3)=90-rad2deg(atan( (mag(1)*cos(tmp_attiNO(1))+mag(3)*sin(tmp_attiNO(1))) /...
       (mag(2)*cos(tmp_attiNO(2))+mag(1)*sin(tmp_attiNO(2))*sin(tmp_attiNO(1))-mag(3)*sin(tmp_attiNO(2))*cos(tmp_attiNO(1)))  ));
    %%  UKF�Ա�
    tic;
    [ukf_q,ukf_X,ukf_p] = USQUE(a,f,Q,R,n,lambda,T,ukf_X,ukf_q,Wibb,[mag;Fb],ukf_p);
    USQUE_t(fix(t*100+1))=toc;
%     [ukf_q,ukf_X,pred_q,pred_mean,pred_cov,aleph] = USQUE2(pred_mean,pred_cov,aleph,pred_q,T,ukf_X,ukf_q,Wibb,[mag;Fb],noise_gyro,noise_wARW,noise_accel,noise_mag);
    ukf_atti=atti_compute(ukf_q(1:4));  
    
    %% ���ݴ洢
    TraceData=[TraceData;t,atti',atti_fromMeasure',atti_fromGyro',ukf_atti'];
    %% ������
    UKFerro=qAntiMatrix([ukf_q(1);-ukf_q(2:4)])*Ture_q;
    ErroData=[ErroData;t,2*acos(abs(UKFerro(1)))*180/pi];
    wibbdata=[wibbdata;t,norm(deg2rad(ukf_X(4:6))-Gyro_r)];
    %%
    t=t+T;
end

RMSE_USQUE=sqrt(sum(ErroData(281:end,2).^2)/((t_stop-2.81)*100))
%% ��̬��ʾ
figure(1);
% unbiased
subplot(3,1,1);plot(TraceData(:,1),TraceData(:,2));hold on;
subplot(3,1,2);plot(TraceData(:,1),TraceData(:,3));hold on;
subplot(3,1,3);plot(TraceData(:,1),TraceData(:,4));hold on;
% atti from measurement vector
subplot(3,1,1);plot(TraceData(:,1),TraceData(:,5));hold on;
subplot(3,1,2);plot(TraceData(:,1),TraceData(:,6));hold on;
subplot(3,1,3);plot(TraceData(:,1),TraceData(:,7));hold on;
% atti from Gyro
subplot(3,1,1);plot(TraceData(:,1),TraceData(:,8));hold on;
subplot(3,1,2);plot(TraceData(:,1),TraceData(:,9));hold on;
subplot(3,1,3);plot(TraceData(:,1),TraceData(:,10));hold on;
% USQUE
subplot(3,1,1);plot(TraceData(:,1),TraceData(:,11));hold off;
subplot(3,1,2);plot(TraceData(:,1),TraceData(:,12));hold off;
subplot(3,1,3);plot(TraceData(:,1),TraceData(:,13));hold off;
legend('True','Gyro','Measurement','USQUE');
%% error
figure(2);
semilogy(ErroData(:,1),ErroData(:,2));hold off;%UKF
legend('USQUE');
figure(3);
semilogy(ErroData(:,1),wibbdata(:,2));
