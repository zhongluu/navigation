% Copyright(c) 2021, by Yulu Zhong. All rights reserved.
% Key Laboratory of Micro-Inertial Instrument and Advanced Navigation Technology of Ministry of Education,
% Southeast University, NanJing, P.R.China 10/31/2021
% based on psins toolbox from http://www.psins.org.cn/
% version:psins210406.rar
% or psins210522.rar
% Acknowledge: Gongmin Yan.
close all;
clear;
glvs;

trj = trjfile('trj10ms.mat');
[nn, ts, nts] = nnts(1, trj.ts);
% sensor's deviation defined -- same as reference
deviationOfGPS = 5; % m
deviationOfGyroV = rad2deg(2.9089e-7 * 60); % rad/sqrt(Hz) --> deg/sqrt(h)
deviationOfGyroU = 9.1989e-7; % rad/s/sqrt(Hz)
deviationOfAccV = 9.8100e-5 / glv.ug; % m/s^2/sqrt(Hz) --> ug/sqrt(Hz)
deviationOfAccU = 6e-5; % m/s^2/s^2/sqrt(Hz)
initialAttiErr = 15 * 60 * ones(1, 3); % deg --> arcmin
initialVelErr = zeros(1, 3);
initialPosErr = zeros(1, 3);
initialBiasG = 10; % deg/hr
initialBiasA = 0.003 / glv.ug; % m/s^2 --> ug

%% init
imuerr = imuerrset(0, 0, deviationOfGyroV, deviationOfAccV); 
imu = imuadderr(trj.imu, imuerr);

% init usque
davp0 = avperrset(0, 0, 0);
% davp0 = avperrset(initialAttiErr, initialVelErr, initialPosErr);
initavp = avpadderr(trj.avp0, davp0);
opt_q = a2qua(initavp(1:3)); % the actual quaternion for filter
pos = initavp(7:9);
vel = initavp(4:6);
esteth = ethinit(pos, vel);
biasG = zeros(3, 1);
biasA = zeros(3, 1);
deltaS = zeros(3, 1); % nominal quaternion
X = [deltaS', vel', pos', biasG', biasA']'; % a 15 state -- aleph0 pos after vel
% prealloc
len = length(imu);
[estRes, estRotationAngle, realRes] = prealloc(fix(len / nn), length(X) + 1, 1, 9);
timebar(nn, len, 'MineAlignUSQUETestNV.');
ki = 1;
for k = 1:nn:len - nn + 1
    k1 = k + nn - 1;
    wvm = imu(k:k1, 1:6); t = imu(k1, end);
%     [phim, dvbm] = cnscl(wvm, 0);

    vn01 = X(4:6); pos01 = X(7:9);
    esteth = ethupdate(esteth, pos01, vn01);
    
    % velocity  updating
    fb = wvm(4:6)' / nts;
    estfn = qmulv(opt_q, fb - X(13:15));
    estan = estfn + esteth.gcc;
    X(4:6) = X(4:6) + estan * nts;
    % position updating
    estMpv = [0, 1/esteth.RMh, 0; 1/esteth.clRNh, 0, 0; 0, 0, 1];
    X(7:9) = X(7:9) + estMpv*X(4:6) * nts;
     % attitude updating
    est_omega = wvm(1:3)' - X(10:12) * nts;
    opt_q = qupdt2(opt_q, est_omega, esteth.wnin * nts);   
    %% other
    realQ = a2qua(trj.avp(k1, 1:3));
    erroQ = qAntiMatrix(realQ) * [opt_q(1); -opt_q(2:4)];   
    % save result
    estAtti = rad2deg(q2att(opt_q)); % deg, m/s, m
    estRes(ki, :) = [estAtti; X(4:end); t]';
    realRes(ki, :) = [rad2deg(trj.avp(k1, 1:3)), trj.avp(k1, 4:9)]; % deg, m/s, m
    estRotationAngle(ki) = 2 * acos(abs(erroQ(1))) * 180 / pi;
    ki = ki + 1;
    timebar;
end
% free unnecessary space for ram
estRes(ki:end,:) = []; realRes(ki:end,:) = []; estRotationAngle(ki:end,:) = []; 
%% plot
% attitude
figure('Name', 'attitude');
subplot(3, 1, 1); plot(estRes(:, end), [realRes(:, 1), estRes(:, 1)]); hold on; % pitch
subplot(3, 1, 2); plot(estRes(:, end), [realRes(:, 2), estRes(:, 2)]); hold on; % yaw
subplot(3, 1, 3); plot(estRes(:, end), [realRes(:, 3), estRes(:, 3)]); hold off; % roll
legend('True','Estimation');
% velocity
figure('Name', 'velocity');
subplot(3, 1, 1); plot(estRes(:, end), [realRes(:, 4), estRes(:, 4)]); hold on; % east
subplot(3, 1, 2); plot(estRes(:, end), [realRes(:, 5), estRes(:, 5)]); hold on; % north
subplot(3, 1, 3); plot(estRes(:, end), [realRes(:, 6), estRes(:, 6)]); hold on; % up
legend('True','Estimation');
% position
figure('Name', 'position');
subplot(3, 1, 1); plot(estRes(:, end), [realRes(:, 7), estRes(:, 7)]); hold on; % latitude
subplot(3, 1, 2); plot(estRes(:, end), [realRes(:, 8), estRes(:, 8)]); hold on; % longitude
subplot(3, 1, 3); plot(estRes(:, end), [realRes(:, 9), estRes(:, 9)]); hold on; % height
legend('True','Estimation');
% rotation angle
figure('Name', 'rotation angle from error quaternion');
semilogy(estRes(:, end), estRotationAngle); hold off; %UKF
legend('error');
% trajectory
figure('Name', 'vehicle tracjectory');
plot3([realRes(:, 7), estRes(:, 7)],[realRes(:, 8), estRes(:, 8)],[realRes(:, 9), estRes(:, 9)]);
%% other needed function
