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
 
clear;
close all;
clc;
addpath('basic functions');
%% ����ʱ������
n=0;
t=0;
T=0.01;
Nsind=0;   
M=500;%������
t_stop=10;
flag=0;
%% ����������
atti=zeros(3,1);     %��������������򣨵�λ���ȣ�
atti_rate=zeros(3,1);%������ʡ��������ʡ��������ʣ���λ����/�룩

atti(1,1)=-0.2420;
atti(2,1)=1.5754;
atti(3,1)=82.9399;      %��ʼ����Ƕȣ���λ���ȣ�
%% Ԥ���估��ʼ����̬����
attiN=atti;              %��ʼ��̬�뺽����̬һ�£������ó�ʼ��׼�����滻��
%% Ԥ����IMU���
Wibb=zeros(3,1);    %����ϵ���������   ����λ����/�룩
Fb=zeros(3,1);      %����ϵ���ٶȼ���� ����λ����/��/�룩
Gyro_r=zeros(3,1);  % ����һ������ɷ���̣�����/�룩
Acc_r =zeros(3,1);  % ���ٶ�һ������ɷ���̣���/��/�룩
real_Wibb=zeros(3,1);
%%  ��ʼ����Ԫ��
q= Eluer_to_quaternion( attiN );
%% �����˲���ز���
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
%% sensors specification ��OpenIMU300ZI��
noise_wARW=deg2rad(50e-4);
noise_gyro=deg2rad(70e-4);% rad/s
noise_accel=10e-6;% g
noise_mag=25e-9;% T

%% sensors specification ��OpenIMU300ZI:mag_z1.csv��
noise_wARW=[deg2rad(61e-4) 0 0;0 deg2rad(51e-4) 0;0 0 deg2rad(55e-4)];
noise_gyro=[deg2rad(73e-4) 0 0;0 deg2rad(70e-4) 0;0 0 deg2rad(67e-4)]';% rad/s
noise_fARW=[(10e-4) 0 0;0 (10e-4) 0;0 0 (10e-4)];
noise_accel=[4.25e-4 0 0;0 5.85e-4 0;0 0 6.38e-4];
noise_mag=2.5e-8*eye(3);% T
% %% sensors specification ��OpenIMU300ZI:static2_z1.csv ת̨�ɼ���
noise_wARW=[deg2rad(67e-4) 0 0;0 deg2rad(50e-4) 0;0 0 deg2rad(54e-4)];
noise_gyro=[deg2rad(41e-4) 0 0;0 deg2rad(79e-4) 0;0 0 deg2rad(70e-4)]'.*5;% rad/s
noise_fARW=[(2.68e-4) 0 0;0 (2.8607e-4) 0;0 0 (2.9843e-4)];
noise_accel=[4e-4 0 0;0 4.4e-4 0;0 0 9.68e-4];
noise_mag=[8.32e-8 0 0;0 9.08e-8 0;0 0 1.45e-7];% T
noise_mag=[3.31537718833425e-08,0,0;0,7.02234382409494e-08,0;0,0,5.87903320244251e-08];
%% ��ʼ��PF
inital_t=0;
PF_X=zeros(4,M);
PF_estX=zeros(4,1);
PF_Weight=ones(1,M)*(1/M);
PF_wibbARW=zeros(3,M);
PF_estWibb=zeros(3,1);
%% �������ݴ洢�ռ�
pqpfparticles=[];
TraceData=[];
RegKGData=[];
wibbdata=[];
real_qData=[];
sq_Data=[];
ErroData=[];
boundary_data=[];
Wibb_gyroscope=[];Acceleration=[];Magnetitude=[];Flag_disturbance=[];
%% ��ʼ����������̬��ز���
nosi_q=q;
nosi_atti = atti;
%% ��ʼ����������̬��ز���
real_q = q;
real_atti = atti;
atti_d = [];

while t<=t_stop
    %% �������ݲ���
    [t,atti,atti_rate] = traceset2(t,T,atti,atti_rate);%�����������������й켣����
    [ Wibb,Gyro_r ] = gyroscope(t,T,atti,atti_rate,Gyro_r,noise_wARW,noise_gyro);%�������������������
    [ real_Wibb] = real_gyroscope(atti,atti_rate);%�������������������
    [Fb,Acc_r] = accelerometer(t,T,atti,noise_accel,Acc_r);%���������ٶȼ��������
    [mag] = magnet(t,atti,noise_mag);
    Wibb_gyroscope = [Wibb_gyroscope,Wibb];
    Acceleration = [Acceleration,Fb];
    Magnetitude = [Magnetitude,mag];
    %% ��������̬���
    if t == 0
        real_q = q;
    else
        [real_q] = quaternion_updata(T,real_Wibb,real_q);% ���������ٶȻ���
    end
    real_atti = atti_compute(real_q);
    %%
    [Ture_q] = Eluer_to_quaternion( atti );
    %% ��������̬���
    [nosi_q] = quaternion_updata(T,Wibb,nosi_q);% ���������ٶȻ���
    nosi_atti = atti_compute(nosi_q);
    [ attiNO ] = acc_cal_atti( Fb );%%�Ƕ�
    tmp_attiNO = deg2rad(attiNO);
    attiNO(3) = 90-rad2deg(atan( (mag(1)*cos(tmp_attiNO(1))+mag(3)*sin(tmp_attiNO(1))) /...
        (mag(2)*cos(tmp_attiNO(2))+mag(1)*sin(tmp_attiNO(2))*sin(tmp_attiNO(1))-mag(3)*sin(tmp_attiNO(2))*cos(tmp_attiNO(1)))  ));
    nosi_q_1 = Eluer_to_quaternion(attiNO);
%     %% PF
    if t==0
        [PF_X] = particle_init_mine(mag,[0.48e-3;0.0;0.48e-3],M);
    end
    if t<=2.80
         [ PF_estX,PF_X,Neff,PF_Weight,PF_wibbARW] = initalparticleswithoutKLD(t,T,PF_X,Wibb,Fb,mag,PF_Weight,PF_wibbARW,noise_wARW,noise_gyro,noise_accel,noise_mag );
         Wibb_last=Wibb;WibbGyro_r_last=Gyro_r;
    else
    [ PF_estX,PF_X,PF_Weight,PF_wibbARW,flag] = SPF_of_atti2(T,M,PF_X,PF_Weight,Wibb_last-rad2deg(WibbGyro_r_last),Wibb,Fb,mag,PF_wibbARW,noise_wARW,noise_gyro,noise_mag,noise_accel);
    Wibb_last = Wibb;WibbGyro_r_last = Gyro_r;
    Flag_disturbance = [Flag_disturbance,flag];
    end
    [PF_atti] = atti_compute(PF_estX);
    PF_estWibb = PF_wibbARW * PF_Weight';
    for m=1:M
        tmp_attiN(:,m) = atti_compute(PF_X(:,m)) - atti;
    end
    upboundary_atti = +3 * std(tmp_attiN,0,2); 
    downboundary_atti = -3 * std(tmp_attiN,0,2);
    boundary_data = [boundary_data;t, upboundary_atti', downboundary_atti'];
    %% ���ݴ洢
    TraceData = [TraceData; t, atti', attiNO', nosi_atti', real_atti', PF_atti'];
    wibbdata = [wibbdata; t, norm(PF_estWibb - Gyro_r)];
    real_qData = [real_qData'; real_q']';
    %% ������
    PFerro = qAntiMatrix([PF_estX(1); - PF_estX(2:4)]) * Ture_q;
    Measurement_erro = qAntiMatrix([nosi_q_1(1); - nosi_q_1(2:4)]) * Ture_q;
    Gyro_erroWithoutNoiseInte = qAntiMatrix([real_q(1); - real_q(2:4)]) * Ture_q; %%���������ٶȻ���  
    Gyro_erro = qAntiMatrix([nosi_q(1); - nosi_q(2:4)]) * Ture_q; %%���������ٶȻ���
    ErroData = [ErroData;t, 2*acos(abs(PFerro(1)))*180 / pi,...
                2 * acos(abs(Measurement_erro(1)))*180 / pi,...
                2 * acos(abs(Gyro_erroWithoutNoiseInte(1))) * 180 / pi,...
                2 * acos(abs(Gyro_erro(1))) * 180 / pi];
    %%
    t=t+T;
end
RMSE_PF = sqrt(sum(ErroData(500:end,2).^2) / ((t_stop-5)*100));
RMSE_Mea = sqrt(sum(ErroData(1:end,3).^2) / ((t_stop)*100));
RMSE_GyroWithoutNoise = sqrt(sum(ErroData(1:end,4).^2) / ((t_stop)*100));
RMSE_Gyro = sqrt(sum(ErroData(1:end,5).^2) / ((t_stop)*100));

% %% ��̬��ʾ
figure(1);
subplot(3,1,1);plot(TraceData(:,1),TraceData(:,2));hold on;%������-�����
subplot(3,1,2);plot(TraceData(:,1),TraceData(:,3));hold on;%������-������
subplot(3,1,3);plot(TraceData(:,1),TraceData(:,4));hold on;%������-�����

subplot(3,1,1);plot(TraceData(:,1),TraceData(:,5));hold on;%���ٶȼ�-�����
subplot(3,1,2);plot(TraceData(:,1),TraceData(:,6));hold on;%���ٶȼ�-������
subplot(3,1,3);plot(TraceData(:,1),TraceData(:,7));hold on;%���ٶȼ�-�����

subplot(3,1,1);plot(TraceData(:,1),TraceData(:,8));hold on;%���������������-�����
subplot(3,1,2);plot(TraceData(:,1),TraceData(:,9));hold on;%���������������-������
subplot(3,1,3);plot(TraceData(:,1),TraceData(:,10));hold on;%���������������-�����

subplot(3,1,1);plot(TraceData(:,1),TraceData(:,11));hold on;%����������ʵ�����������-�����
subplot(3,1,2);plot(TraceData(:,1),TraceData(:,12));hold on;%����������ʵ�����������-������
subplot(3,1,3);plot(TraceData(:,1),TraceData(:,13));hold on;%����������ʵ�����������-�����

subplot(3,1,1);plot(TraceData(:,1),TraceData(:,14));hold off;%PF-�����
subplot(3,1,2);plot(TraceData(:,1),TraceData(:,15));hold off;%PF-������
subplot(3,1,3);plot(TraceData(:,1),TraceData(:,16));hold off;%PF-�����

legend('True','Fb','Gyro','Real','PF');
%% error
figure(2);
semilogy(ErroData(:,1),ErroData(:,2));hold on;
semilogy(ErroData(:,1),ErroData(:,3));hold on;
semilogy(ErroData(:,1),ErroData(:,4));hold on;
semilogy(ErroData(:,1),ErroData(:,5));hold off;
legend('QPF','Measurement','Gyro_WithoutNoise','Gyro_Noise');
figure(3);
subplot(3,1,1);plot(TraceData(:,1),TraceData(:,14)-TraceData(:,2),'b');hold on;
subplot(3,1,1);plot(boundary_data(:,1),boundary_data(:,2),'r',boundary_data(:,1),boundary_data(:,5),'r');hold off;axis([0,t,-0.05,0.05]);
ylabel('Roll(Deg)');
subplot(3,1,2);plot(TraceData(:,1),TraceData(:,15)-TraceData(:,3),'b');hold on;
subplot(3,1,2);plot(boundary_data(:,1),boundary_data(:,3),'r',boundary_data(:,1),boundary_data(:,6),'r');hold off;axis([0,t,-0.05,0.05]);
ylabel('Pitch(Deg)');
subplot(3,1,3);plot(TraceData(:,1),TraceData(:,16)-TraceData(:,4),'b');hold on;
subplot(3,1,3);plot(boundary_data(:,1),boundary_data(:,4),'r',boundary_data(:,1),boundary_data(:,7),'r');hold off;axis([0,t,-0.05,0.05]);
ylabel('Yaw(Deg)');
xlabel('Time(Sec)');
figure(4);
semilogy(ErroData(:,1),wibbdata(:,2));