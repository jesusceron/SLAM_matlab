function [positions, beac_rssi_fixed, beac_rssi_fixed_filtered, beac_rssi_activity, beac_rssi_activity_filtered] = Load_data()

%% Load files
Acc_Gyr_Beac_data= readtable('E:\Google Drive\IL_HAR_app_Matlab\dataset_synchronized\11.csv');
fs = 204.8; % IMU sample rate
stance_phase = load('E:\Google Drive\IL_HAR_app_Matlab\dataset_synchronized\stance_phase_p11.mat');


%% THIS CHANGES FOR EACH PARTICIPANT
gyro_i_c_s = 1; % initial calibration sample
samples_c = 3000;  % asumo 5 segundos parado (y fs=204.8 Hz)
start_sample = 3900;
final_sample = 62720;

%% Load data.

% Load IMU data
acc_complete = table2array(Acc_Gyr_Beac_data(:,23:25));
gyro_complete = table2array(Acc_Gyr_Beac_data(:,26:28));

% Load BEACONS data.
beac_rssi_fixed = table2array(Acc_Gyr_Beac_data(:,...
    {'RSSIs_beacon_1',...   % Coconut:  Room
    'RSSIs_beacon_2',...    % Mint:     Kitchen
    'RSSIs_beacon_3',...    % Ice:      Bathroom
    'RSSIs_beacon_4',...    % Blueber:  Dining room
    'RSSIs_beacon_5'}));    % P2:       Living room

beac_rssi_activity = table2array(Acc_Gyr_Beac_data(:,...
    {'RSSIs_beacon_6',...   % P1:       Door
    'RSSIs_beacon_7',...    % B2:       Toilet lid
    'RSSIs_beacon_8',...    % B1:       Broom
    'RSSIs_beacon_9',...    % G2:       Pitcher
    'RSSIs_beacon_10'}));   % G1:       Hair brush

beac_motion = table2array(Acc_Gyr_Beac_data(:,{'motion_state_beacon_6',...
    'motion_state_beacon_7',...
    'motion_state_beacon_8',...
    'motion_state_beacon_9',...
    'motion_state_beacon_10'}));

%% Data pre-processing
w = gyro_i_c_s : gyro_i_c_s + samples_c;

acc = zeros(length(acc_complete),3);
acc(:,1) = acc_complete(:,2); 
acc(:,2) = acc_complete(:,3); % y = z
acc(:,3) = acc_complete(:,1); 

gyr = zeros(length(gyro_complete),3);
gyr(:,1) = gyro_complete(:,2); 
gyr(:,2) = gyro_complete(:,3); 
gyr(:,3) = gyro_complete(:,1); 
gyr = deg2rad(gyr);

gravity=mean(sqrt(acc(w,1).^2 + acc(w,2).^2 + acc(w,3).^2));  % m/s^2
acc_mean=[mean(acc(w,1)),mean(acc(w,2)),mean(acc(w,3))];

% Remove bias Gyro
bias_gyr=[mean(gyr(w, 1)), mean(gyr(w, 2)), mean(gyr(w, 3))];
gyr_unbiased=[gyr(:,1) - bias_gyr(1),...
    gyr(:,2) - bias_gyr(2),...
    gyr(:,3) - bias_gyr(3)];

acc = acc(start_sample:final_sample,:);
gyr_unbiased = gyr_unbiased(start_sample:final_sample,:);
beac_rssi_fixed = beac_rssi_fixed(start_sample:final_sample,:);
beac_rssi_activity = beac_rssi_activity(start_sample:final_sample,:);
beac_motion = beac_motion(start_sample:final_sample,:);
%beac_rssi_fixed = [zeros(30,5); beac_rssi_fixed(start_sample:final_sample - 30,:)];
%beac_rssi_activity = [zeros(30,5); beac_rssi_activity(start_sample:final_sample - 30,:)];
%beac_motion = beac_motion(start_sample:final_sample - 30,:);

% Filtering the beacons data in Python
[beac_rssi_fixed_filtered, beac_rssi_activity_filtered] = rssiKF(beac_rssi_fixed,beac_rssi_activity);


%% -----Step detection---------------------------
idx_fig=20;
% [Num_steps,Step_events,StancePhase,idx_fig]=StepDetection_Acel(acc,1,idx_fig);

Step_events = stance_phase.stance_phase(:,1)';


%% -------------------- Trajectory reconstruction ZUPT KF -----------------------%
[positions]=ZUPT_KF(acc, gyr_unbiased, gravity, acc_mean, fs, Step_events, idx_fig+4);


    