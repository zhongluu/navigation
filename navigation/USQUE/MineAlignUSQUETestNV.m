% Copyright(c) 2021, by Yulu Zhong. All rights reserved.
% Key Laboratory of Micro-Inertial Instrument and Advanced Navigation Technology of Ministry of Education,
% Southeast University, NanJing, P.R.China 10/31/2021
% based on psins toolbox from http://www.psins.org.cn/
% version:psins210406.rar
% or psins210522.rar
% Acknowledge: Gongmin Yan and Kailong Li.
close all;
clear;
glvs;
psinstypedef(153);
trj = trjfile('trj10ms.mat');
[nn, ts, nts] = nnts(1, trj.ts);
nts2 = nts / 2;
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
imuerr = imuerrset(initialBiasG, initialBiasA, deviationOfGyroV, deviationOfAccV); % same as reference
imu = imuadderr(trj.imu, imuerr);
%% a direct sigma point filter : Kg and Ka are omitted
% reference: JOHN L. CRASSIDIS. "Sigma-point_Kalman_filtering_for_integrated_GPS_and_inertial_navigation"
%            Crassidis, John L, Markley, F. "Landis. Unscented Filtering for Spacecraft Attitude Estimation"
%            Rudolph, Vdm, Wan, E., Julier, S. "Sigma-Point Kalman Filters for Nonlinear Estimation and Sensor Fusion: Applications to Integrated Navigation"
%            John L. Crassidis. "Sigma-Point Kalman Filtering for Integrated GPS and Inertial Navigation"
tmpan = zeros(3, 1); tmpMpv = zeros(3, 3);
% init usque
davp0 = avperrset(initialAttiErr, initialVelErr, initialPosErr);
initavp = avpadderr(trj.avp0, davp0);
opt_q = a2qua(initavp(1:3)); % the actual quaternion for filter
pos = initavp(7:9);
vel = initavp(4:6);
esteth = ethinit(pos, vel);
biasG = zeros(3, 1);
biasA = zeros(3, 1);
deltaS = zeros(3, 1); % nominal quaternion
X = [deltaS', vel', pos', biasG', biasA']'; % a 15 state -- aleph0 pos after vel
% filter required parameter
a = 1;
f = 2 * (a + 1);
% 
P = diag([(deg2rad(15/3))^2 * ones(1, 3), (200/3)^2 * ones(1, 2), (10/3)^2, 1e-6^2 * ones(1, 2), (20/3)^2, ...
        (initialBiasG * glv.dph)^2 * ones(1, 3), (0.005/3)^2 * ones(1, 3)]); % same as reference
% Q0 =  [ (deg2rad(deviationOfGyroV) / 60)^2 * eye(3), zeros(3,9); ...
%      zeros(3, 3), (deviationOfAccV * glv.ug)^2 * eye(3), zeros(3, 6); ...
%                         zeros(6, 12)                    ]; % same as reference
G = [eye(3),zeros(3,9);
    zeros(3),eye(3),zeros(3,6);
    zeros(9,12)];
% Q =  [ (deg2rad(deviationOfGyroV) / 60)^2 * eye(3), zeros(3,12); ...
%      zeros(3, 3), (deviationOfAccV * glv.ug)^2 * eye(3), zeros(3, 9); ...
%                         zeros(9, 15)                    ]; % 
Q = zeros(size(P));
R = diag([(deviationOfGPS/glv.Re)^2*ones(1,2),deviationOfGPS^2]); % same as reference
n = length(X);
kappa = 4 - n; % same as reference
alpha = 1; % same as reference
lambda = alpha^2 * (n + kappa) - n; % same as reference
beta = 2; % same as reference
% prealloc
len = length(imu);
[estRes, estRotationAngle, realRes] = prealloc(fix(len / nn), length(X) + 1, 1, length(X));
[sigmaLine3] = prealloc(fix(len / nn),2*length(X));
timebar(nn, len, 'MineAlignUSQUETestNV.');
ki = 1;
for k = 1:nn:len - nn + 1
    k1 = k + nn - 1;
    wvm = imu(k:k1, 1:6); t = imu(k1, end);
    fb = wvm(4:6)' / nts;

    %% start
    % earth updating
    vn01 = X(4:6); pos01 = X(7:9);
    esteth = ethupdate(esteth, pos01, vn01);
% TODO: some problems for Q    
%     ins.eth = esteth; 
%     ins.nts = nts;
%     ins.vn = X(4:6);
%     ins.Cnb = q2mat(opt_q);
%     [ins.tauG, ins.tauA] = setvals(inf(3,1));
%     ins.fn = qmulv(opt_q, fb);
%     ins.Mpv = [0, 1/esteth.RMh, 0; 1/esteth.clRNh, 0, 0; 0, 0, 1];
%     fhi = kffk(ins);
%     Q =eye(15,15) + G * Q0 * G' * nts;
%     Q = sylvester(fhi,fhi,fhi*Q);
    %% USQUE start
    sigma = [zeros(size(X)), -1 * chol((n + lambda) * (P + Q))', chol((n + lambda) * (P + Q))'];
    sigmadots = length(sigma);
    aleph = repmat(X, 1, sigmadots) + sigma;

    % convert MRPs to error quaternion
    delta_q = zeros(4, sigmadots);
    for i = 1:sigmadots
        scalar_delta_q = (-a * norm(aleph(1:3, i))^2  + f * sqrt(f^2 + (1 - a^2) * norm(aleph(1:3, i))^2)) / (f^2 + norm(aleph(1:3, i))^2);
        vector_delta_q = (1 / f) * (a + scalar_delta_q) .* aleph(1:3, i);
        delta_q(:, i) = [scalar_delta_q; vector_delta_q];
    end

    %% prediction
    sq = zeros(4, sigmadots);
    for i = 1:sigmadots
        sq(:, i) = qmul(delta_q(:, i), opt_q);
    end
     
    % velocity and position updating
    for i = 1:sigmadots
        % velocity
        estfn = qmulv(sq(:, i), (fb - aleph(13:15, i) - (deviationOfAccV * glv.ug) .* randn(3,1)));%
        estan = estfn + esteth.gcc;
        aleph(4:6, i) = aleph(4:6,i) + estan * nts;
        % position
        estMpv = [0, 1/esteth.RMh, 0; 1/esteth.clRNh, 0, 0; 0, 0, 1];
        aleph(7:9, i) = aleph(7:9, i) + estMpv * aleph(4:6,i) * nts;
    end
 
    % attitude updating
    est_omega = repmat(wvm(1:3)', 1, sigmadots) - aleph(10:12, :) * nts - (deg2rad(deviationOfGyroV) / 60) .* randn(3,sigmadots) * nts;% 
    pred_q = zeros(4, sigmadots);
    for i = 1:sigmadots
        pred_q(:, i) = qupdt2(sq(:, i), est_omega(:, i), esteth.wnin * nts);
    end 
   
    pred_delta_q = zeros(4, sigmadots);
    for i = 1:sigmadots
        pred_delta_q(:, i) = qmul(pred_q(:, i), [pred_q(1, 1); -pred_q(2:4, 1)]);
    end
    % bias updating
        
    %convert error quaternion to MRPs
    % aleph(1:3, 1) = zeros(3,1);
    for i = 1:sigmadots
        aleph(1:3, i) = f * (pred_delta_q(2:4, i) ./ (a + pred_delta_q(1, i)));
    end

    % standard SPKF -- deriv Pxx
    pred_mean = (1 / (n + lambda)) * (lambda * aleph(:, 1) + (1/2) * sum(aleph(:, 2:end), 2));
    pred_cov = (1 / (n + lambda)) * (lambda * (aleph(:, 1) - pred_mean) * (aleph(:, 1) - pred_mean)' ...
        + (1/2) * ((aleph(:, 2:end) - repmat(pred_mean, 1, sigmadots - 1)) ...
        * (aleph(:, 2:end) - repmat(pred_mean, 1, sigmadots - 1))')) ...
        + (1-alpha^2 + beta)*((aleph(:, 1) - pred_mean) * (aleph(:, 1) - pred_mean)') + Q;
    %% correction
    if mod(t,1)==0 
        % GPS pos simulation with some white noise
        posGPS = trj.avp(k1, 7:9)' + diag([deviationOfGPS / glv.Re * ones(2, 1); deviationOfGPS]) * randn(3, 1);
        H = [zeros(3,6), eye(3), zeros(3,6)];
        gamma = H * aleph;
        % standard SPKF -- deriv Pyy
        obs_mean = (1 / (n + lambda)) * (lambda * gamma(:, 1) + (1/2) * sum(gamma(:, 2:end), 2));
        obs_cov=(1/(n+lambda))*(lambda*(gamma(:,1)-obs_mean)*(gamma(:,1)-obs_mean)'...
            +(1/2)*((gamma(:,2:end)-repmat(obs_mean,1,sigmadots-1))*(gamma(:,2:end)-repmat(obs_mean,1,sigmadots-1))')) ...
            + (1-alpha^2 + beta)*((gamma(:,1)-obs_mean)*(gamma(:,1)-obs_mean)');
%         cross_corr_cov = pred_cov * H';
%         obs_cov = H * pred_cov * H';
        innov_cov = obs_cov + R;
        % standard SPKF -- deriv Pxy
        cross_corr_cov=(1/(n+lambda))*(lambda*(aleph(:,1)-pred_mean)*(gamma(:,1)-obs_mean)'...
        +(1/2)*((aleph(:,2:end)-repmat(pred_mean,1,sigmadots-1))*(gamma(:,2:end)-repmat(obs_mean,1,sigmadots-1))')) ...
        + (1-alpha^2 + beta)*((aleph(:,1)-pred_mean)*(gamma(:,1)-obs_mean)');
        % standard SPKF -- deriv K
        KGain = cross_corr_cov / (innov_cov);
        % measurement updating
        X = pred_mean + KGain * (posGPS - obs_mean);
        P = pred_cov - KGain * innov_cov * KGain';
        
        % estimate the optimal quaternion
        scalar_opt_delta_q = (-a * norm(X(1:3))^2 + f * sqrt(f^2 + (1 - a^2) * norm(X(1:3))^2)) / (f^2 + norm(X(1:3))^2);
        vector_delta_q = (1 / f) * (a + scalar_opt_delta_q) .* X(1:3);
        opt_delta_q = [scalar_opt_delta_q; vector_delta_q];
        opt_q = qmul(opt_delta_q, pred_q(:,1));
        % reset the MRPs
        X(1:3) = zeros(3, 1);
        % measurement updating completing
        % calcuate rotation angle
        realQ = a2qua(trj.avp(k1, 1:3));
        erroQ = qmul(realQ, [opt_q(1); -opt_q(2:4)]);
        % save filter result
        estAtti = rad2deg(q2att(opt_q)); % deg, m/s, m
        estRes(ki, :) = [estAtti; X(4:end); t]';
        sigmaLine3(ki, :) = [(3 * chol(P)*ones(length(X),1))', (-3 * chol(P)*ones(length(X),1))'];
        realRes(ki, :) = [rad2deg(trj.avp(k1, 1:3)), trj.avp(k1, 4:9), imuerr.eb', imuerr.db']; % deg, m/s, m , rad/s, m/s^2
        estRotationAngle(ki) = 2 * acos(abs(erroQ(1))) * 180 / pi;
        ki = ki + 1;
    else
        % only time updating
        X = pred_mean;
        P = pred_cov;
%         P = (P + P')/2;
        scalar_opt_delta_q = (-a * norm(X(1:3))^2 + f * sqrt(f^2 + (1 - a^2) * norm(X(1:3))^2)) / (f^2 + norm(X(1:3))^2);
        vector_delta_q = (1 / f) * (a + scalar_opt_delta_q) .* X(1:3);
        opt_delta_q = [scalar_opt_delta_q; vector_delta_q];
        opt_q = qmul(opt_delta_q, pred_q(:,1));
%         X(1:3) = zeros(3, 1);
    end
    
    timebar;
end
% free unnecessary space for ram
estRes(ki:end,:) = []; realRes(ki:end,:) = []; estRotationAngle(ki:end,:) = []; sigmaLine3(ki:end,:) = [];
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
% gyro bias estimation
figure('Name', 'gyroscope bias');
subplot(3, 1, 1); plot(estRes(:, end), [realRes(:, 10), estRes(:, 10)]); hold on; % x
subplot(3, 1, 2); plot(estRes(:, end), [realRes(:, 11), estRes(:, 11)]); hold on; % y
subplot(3, 1, 3); plot(estRes(:, end), [realRes(:, 12), estRes(:, 12)]); hold on; % z
legend('True','Estimation');
% acc bias estimation
figure('Name', 'accelerate bias');
subplot(3, 1, 1); plot(estRes(:, end), [realRes(:, 13), estRes(:, 13)]); hold on; % x
subplot(3, 1, 2); plot(estRes(:, end), [realRes(:, 14), estRes(:, 14)]); hold on; % y
subplot(3, 1, 3); plot(estRes(:, end), [realRes(:, 15), estRes(:, 15)]); hold on; % z
legend('True','Estimation');
% rotation angle
figure('Name', 'rotation angle from error quaternion');
semilogy(estRes(:, end), estRotationAngle); hold off; %UKF
legend('USQUE');
% trajectory
figure('Name', 'vehicle tracjectory');
plot3([realRes(:, 7), estRes(:, 7)],[realRes(:, 8), estRes(:, 8)],[realRes(:, 9), estRes(:, 9)]);
legend('True','Estimation');