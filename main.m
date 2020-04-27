clc;
clearvars -except positions step_events beac_rssi_fixed_filtered beac_rssi_activity_filtered beac_motion
close all;

participant = 1;

N_PARTICLES = 500;
N_LM = 10;
LM_SIZE = 2;
THRESHOLD_RESAMPLE = N_PARTICLES / 2;
STD_PERSON_POSITION = 0.1;

% [positions, step_events,beac_rssi_fixed_filtered, beac_rssi_activity_filtered, beac_motion] = Load_data(participant);
% step_events = [1, step_events, length(beac_motion)]; % Add '1' and the lenght of the dataset as points of support


rssi_room = beac_rssi_fixed_filtered(:,1);
rssi_kitchen = beac_rssi_fixed_filtered(:,2);
rssi_bathroom = beac_rssi_fixed_filtered(:,3);
rssi_dining = beac_rssi_fixed_filtered(:,4);
rssi_living = beac_rssi_fixed_filtered(:,5);

rssi_door = beac_rssi_activity_filtered(:,1);
rssi_toilet = beac_rssi_activity_filtered(:,2);
rssi_broom = beac_rssi_activity_filtered(:,3);
rssi_pitcher = beac_rssi_activity_filtered(:,4);
rssi_brush = beac_rssi_activity_filtered(:,5);

door_moving = beac_motion(:,1);
toilet_moving = beac_motion(:,2);
broom_moving = beac_motion(:,3);
moving_window = 0;
pitcher_moving = beac_motion(:,4);
brush_moving = beac_motion(:,5);


%% CREATION OF PARTICLES
particles = [];
trajectories = cell(1,N_PARTICLES);
particles_weights = ones(1,N_PARTICLES) / N_PARTICLES;


figure
hold on
for i_particle=1:N_PARTICLES
    particles = [particles, Particle(N_PARTICLES,N_LM,LM_SIZE)];
    trajectory = plot(nan,nan);
    trajectories{i_particle} = trajectory;
end
particles = init_beacons_position(particles,N_PARTICLES);



%% MAIN LOOP
x_prev = positions(1, 1);
y_prev = positions(1, 2);
for i_stride=2:length(positions)
    disp(i_stride);
    
    % PREDICTION
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = positions(i_stride, 1);
    y = positions(i_stride, 2);
    for i_particle=1:N_PARTICLES
        
        dist_x = (x - x_prev) + randn() * STD_PERSON_POSITION; % randn()=-1.2849
        dist_y = (y - y_prev) + randn() * STD_PERSON_POSITION; % randn()=-0.1054
        particles(1,i_particle).X = particles(i_particle).X + dist_x;
        particles(1,i_particle).Y = particles(i_particle).Y + dist_y;
        
        particles(1,i_particle).T_x = [particles(1,i_particle).T_x, particles(1,i_particle).X];
        particles(1,i_particle).T_y = [particles(1,i_particle).T_y, particles(1,i_particle).Y];
        
        trajectory_x_data = [particles(1,i_particle).T_x];
        trajectory_y_data = [particles(1,i_particle).T_y];
        set(trajectories{i_particle},'XData',trajectory_x_data, 'YData',trajectory_y_data);
    end
    x_prev = x;
    y_prev = y;
    pause(0.05)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % UPDATE
    initial_sample = step_events(i_stride-1);
    final_sample = step_events(i_stride);
    
    % Do not include RSSI measurements while person is still
    if (final_sample - initial_sample) > 615
        initial_sample = final_sample - 410;
    end
    
    beac_moving = [door_moving(initial_sample:final_sample), toilet_moving(initial_sample:final_sample), broom_moving(initial_sample:final_sample),...
        pitcher_moving(initial_sample:final_sample), brush_moving(initial_sample:final_sample)];
    
    rssis = [rssi_door(initial_sample:final_sample), rssi_toilet(initial_sample:final_sample), rssi_broom(initial_sample:final_sample),...
        rssi_pitcher(initial_sample:final_sample), rssi_brush(initial_sample:final_sample),...
        rssi_room(initial_sample:final_sample), rssi_kitchen(initial_sample:final_sample), rssi_bathroom(initial_sample:final_sample),...
        rssi_dining(initial_sample:final_sample), rssi_living(initial_sample:final_sample)];
    
    for i_rssis=1:size(rssis,1)
        non_zero_indexes = find(rssis(i_rssis,:));
        if ~isempty(non_zero_indexes)
            
            for j_rssi=non_zero_indexes
                rssi = rssis(i_rssis,j_rssi);
                current_sample = (initial_sample -1 +i_rssis);
                
                if rssi > -85
                    
                    switch j_rssi
                        
                        case 3 % broom
                            
                            moving = beac_moving(i_rssis,j_rssi);
                            % Do not update the beacon landmark during 3 second after a movement being detected.
                            if moving
                                moving_window = initial_sample - 1 + i_rssis + 615;
                                continue
                            end
                            if moving_window > current_sample
                                continue
                            end
                            
                            % Beacon is still, update landmark
                            [particles, particles_weights, trajectories] = ...
                                update_resample(particles, rssi, j_rssi, particles_weights, N_PARTICLES, trajectories, THRESHOLD_RESAMPLE);
                            
                        otherwise
                            
                            [particles, particles_weights, trajectories] = ...
                                update_resample(particles, rssi, j_rssi, particles_weights, N_PARTICLES, trajectories, THRESHOLD_RESAMPLE);
                            
                    end
                end
            end
        end
    end
    
end


for i_lm=1:N_LM
    hold on
    plot(particles(1,1).Lm(i_lm,1),particles(1,1).Lm(i_lm,2),'d');
end

function [particles, particles_weights, trajectories] = update_resample(particles, rssi, j_rssi, particles_weights, N_PARTICLES, trajectories, THRESHOLD_RESAMPLE)
    
    % To update
    [particles, particles_weights] = update(particles, rssi, j_rssi, particles_weights, N_PARTICLES);

    % To resample
    neff = 1 / sum(sqrt(particles_weights));
    if neff < THRESHOLD_RESAMPLE
        %indexes = double(py.mifunc.stratified_resample(particles_weights))+1;
        indexes = double(py.mifunc.systematic_resample(particles_weights))+1;
        [particles, particles_weights, trajectories] = resample_from_index(particles, particles_weights, indexes, N_PARTICLES, trajectories);
        %particles_weights = ones(1,N_PARTICLES) / N_PARTICLES; % PARTICLES WEIGHTS ARE INITILIZED TO AND EQUAL VALUE

    end
end

function [particles, particles_weights] = update(particles, rssi, lm_id, particles_weights, N_PARTICLES)
    
    distance_rssi = 10 ^ ((rssi + 83) / -10 * 0.6);  % Calculate distance to beacon
    R = distance_rssi / 10;

    for i_particle=1:N_PARTICLES
        % Get the distance between the person and the landmark
        dx = particles(1,i_particle).X - particles(1,i_particle).Lm(lm_id,1);
        dy = particles(i_particle).Y - particles(i_particle).Lm(lm_id,2);
        dist = sqrt((dx * dx) + (dy * dy));
        residual = distance_rssi - dist;

        % Compute Jacobians
        H = [-dx / dist, -dy / dist];

        % Compute covariance of the residual
        % covV = H * Cov_s * H^T + error
        HxCov = [particles(1,i_particle).LmP{1,lm_id}(1, 1) * H(1) + particles(1,i_particle).LmP{1,lm_id}(1, 2) * H(2),...
                 particles(1,i_particle).LmP{1,lm_id}(2, 1) * H(1) + particles(1,i_particle).LmP{1,lm_id}(2, 2) * H(2)];

        covV = (HxCov(1) * H(1)) + (HxCov(2) * H(2)) + R;

        % Calculate Kalman gain
        K_gain = [HxCov(1) * (1 / covV), HxCov(2) * (1.0 / covV)];

        % Calculate the new landmark position
        lm_x = particles(1,i_particle).Lm(lm_id,1) + (K_gain(1) * residual);
        lm_y = particles(1,i_particle).Lm(lm_id,2) + (K_gain(2) * residual);

        % Calculate the new covariance matrix of the landmark
        % cov_t = cov_t-1 - K * covV * K^T
        lm_P_aux = [[K_gain(1) * K_gain(1) * covV, K_gain(1) * K_gain(2) * covV];...
                    [K_gain(2) * K_gain(1) * covV, K_gain(2) * K_gain(2) * covV]];

        lm_P = [[particles(1,i_particle).LmP{1,lm_id}(1, 1) - lm_P_aux(1,1),...
                 particles(1,i_particle).LmP{1,lm_id}(1, 2) - lm_P_aux(1,2)];...
                [particles(1,i_particle).LmP{1,lm_id}(2, 1) - lm_P_aux(2,1),...
                 particles(1,i_particle).LmP{1,lm_id}(2, 2) - lm_P_aux(2,2)]];

        % Update landmark in particle
        particles(1,i_particle).Lm(lm_id, 1) = lm_x;
        particles(1,i_particle).Lm(lm_id, 2) = lm_y;
        particles(1,i_particle).LmP{1,lm_id} = lm_P;

        % Update particles weight
        particles_weights(i_particle) = py.mifunc.calculate_weight(particles_weights(i_particle), dist, covV, distance_rssi);

        %particles_weights(i_particle) = (particles_weights(i_particle) * ...
        %                                (1 / (covV * sqrt(2*pi))) * exp((-(distance_rssi - dist)^2) / (2 * covV^2))) + 1*exp(-300);
    end
    particles_weights = particles_weights / sum(particles_weights);
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).W = particles_weights(i_particle);
    end
end

function [new_particles, new_particles_weights, trajectories] = resample_from_index(particles, particles_weights, indexes, N_PARTICLES, trajectories)
    new_particles_weights = zeros(1,N_PARTICLES);
    new_particles = particles;
    for i_particle = 1:N_PARTICLES
        new_particles_weights(i_particle) = particles_weights(indexes(i_particle));
        new_particles(1,i_particle) = particles(1,indexes(i_particle));
    end
    
    new_particles_weights = new_particles_weights / sum(new_particles_weights);
    for i_particle=1:N_PARTICLES
        %new_particles(1,i_particle).W =
        %new_particles_weights(i_particle);%%        
        new_particles(1,i_particle).W = 1 / N_PARTICLES;
        new_particles_weights(i_particle) = 1 / N_PARTICLES;%%
        
        trajectory_x_data = [new_particles(1,i_particle).T_x];
        trajectory_y_data = [new_particles(1,i_particle).T_y];
        set(trajectories{i_particle},'XData',trajectory_x_data, 'YData',trajectory_y_data);
        plot(particles(1,1).Lm(i_lm,1),particles(1,1).Lm(i_lm,2),'d');
    end
    pause(0.05)
end











    