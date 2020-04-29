clc;
clearvars -except positions step_events beac_rssi_fixed_filtered beac_rssi_activity_filtered beac_motion trajectory_parts
close all;

participant = 2;
% do not use 4,6,7
N_PARTICLES = 5;
N_LM = 10;
LM_SIZE = 2;
THRESHOLD_RESAMPLE = N_PARTICLES / 2;
STD_PERSON_POSITION = 0.1;


[positions, step_events,beac_rssi_fixed_filtered, beac_rssi_activity_filtered, beac_motion] = Load_data(participant);
step_events = [1, step_events, length(beac_motion)]; % Add '1' and the lenght of the dataset as points of support


% Get strides where the beacons are moved. That indicates that is the
% begining of a new trajectory
[trajectory_parts] = GetTrajectoryParts(beac_motion,step_events);

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


figure
xlim([0,20])
ylim([0, 11])
map
hold on
grid on


%% MAIN LOOP

for i_trajectory=1:size(trajectory_parts,1)
        
    x_initial = trajectory_parts(i_trajectory,3);
    y_initial = trajectory_parts(i_trajectory,4);
    initial_step = trajectory_parts(i_trajectory,1);
    final_step = trajectory_parts(i_trajectory,2);
    
    pos_trajectory_part = positions(initial_step: final_step,:);
        
    correction_x = pos_trajectory_part(1, 1) - x_initial;
    correction_y = pos_trajectory_part(1, 2) - y_initial;
    
    pos_trajectory_part(:, 1:2) = pos_trajectory_part(: ,1:2) - [correction_x correction_y];

    
    % CREATION OF PARTICLES
    particles = [];
    trajectories = cell(1,N_PARTICLES);
    landmarks_figures = cell(1,N_LM);
    particles_weights = ones(1,N_PARTICLES) / N_PARTICLES;
    for i_particle=1:N_PARTICLES
        particles = [particles, Particle(N_PARTICLES,N_LM,LM_SIZE, x_initial, y_initial)];
        trajectory = plot(nan, nan, 'color', [.9 .9 .9]);
        trajectories{i_particle} = trajectory;
    end
    for i_landmark=1:N_LM
        landmark_figure = plot(nan, nan,'d');
        landmarks_figures{i_landmark} = landmark_figure;
    end
    
    particles = init_beacons_position(particles,N_PARTICLES);
    current_best_particle = N_PARTICLES;
    previous_best_particle = N_PARTICLES;
    
    
    
    for i_step=2:size(pos_trajectory_part,1)
        disp(i_step)
        % PREDICTION
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        x = pos_trajectory_part(i_step, 1);
        y = pos_trajectory_part(i_step, 2);
        for i_particle=1:N_PARTICLES
            
            dist_x = (x - x_initial) + randn() * STD_PERSON_POSITION;
            dist_y = (y - y_initial) + randn() * STD_PERSON_POSITION;
            particles(1,i_particle).X = particles(i_particle).X + dist_x;
            particles(1,i_particle).Y = particles(i_particle).Y + dist_y;
            
            particles(1,i_particle).T_x = [particles(1,i_particle).T_x, particles(1,i_particle).X];
            particles(1,i_particle).T_y = [particles(1,i_particle).T_y, particles(1,i_particle).Y];
            
            trajectory_x_data = [particles(1,i_particle).T_x];
            trajectory_y_data = [particles(1,i_particle).T_y];
            
            set(trajectories{i_particle},'XData',trajectory_x_data, 'YData',trajectory_y_data, 'color', [.9 .9 .9]);
            
        end
        
        set(trajectories{current_best_particle},'color', 'b');
        
        % Plot the landmarks of the best particle.
        for i_lm=1:N_LM
            hold on
            set(landmarks_figures{i_lm},'XData',particles(1,current_best_particle).Lm(i_lm,1), 'YData',particles(1,current_best_particle).Lm(i_lm,2));
        end
        
        x_initial = x;
        y_initial = y;
        pause(0.05)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % UPDATE
        initial_sample = step_events(i_step-1);
        final_sample = step_events(i_step);
        % Do not include RSSI measurements while person is still
        if (final_sample - initial_sample) > 1024
            initial_sample = final_sample - 1024;
        end
        
        rssis = [rssi_room(initial_sample:final_sample), rssi_kitchen(initial_sample:final_sample), rssi_bathroom(initial_sample:final_sample),...
            rssi_dining(initial_sample:final_sample), rssi_living(initial_sample:final_sample)];
        
        for i_rssis=1:size(rssis,1)
            non_zero_indexes = find(rssis(i_rssis,:));
            if ~isempty(non_zero_indexes)
                
                for j_rssi=non_zero_indexes
                    rssi = rssis(i_rssis,j_rssi);
                    current_sample = (initial_sample -1 +i_rssis);
                    
                    if rssi >-85
                        [particles, particles_weights, trajectories, current_best_particle] = ...
                            update_resample(particles, rssi, j_rssi+5, particles_weights, N_PARTICLES, trajectories, THRESHOLD_RESAMPLE, previous_best_particle);
                    end
                    
                    
                end
            end
        end
        
        
        
        %             %%%%%% DETECTION OF BEACONS' MOVEMENT TO RESET LOCATION %%%%%%%%%%%%%%
        %             beac_moving = [door_moving(initial_sample:final_sample), toilet_moving(initial_sample:final_sample), broom_moving(initial_sample:final_sample),...
        %                 pitcher_moving(initial_sample:final_sample), brush_moving(initial_sample:final_sample)];
        %
        %             for i_beac_moving=1:size(beac_moving,2)
        %                 if ~isempty(find(beac_moving(:,i_beac_moving), 1))
        %                     switch i_beac_moving
        %                         case 1 % door
        %                             if ~door_already_detected
        %                                 initial_stride = i_stride;
        %                                 correction_x_door = positions(initial_stride, 1) - 5.4;
        %                                 correction_y_door = positions(initial_stride, 2) - 5.5;
        %                                 positions(initial_stride:end,1) = positions(initial_stride:end,1) -correction_x_door;
        %                                 positions(initial_stride:end,2) = positions(initial_stride:end,2) -correction_y_door;
        %                                 x_prev = 5.4;
        %                                 y_prev = 5.5;
        %                                 no_movement = false;
        %                                 door_already_detected = true;
        %                                 break
        %                             end
        %
        %                         case 2 % toilet
        %
        %                             if ~toilet_already_detected
        %                                 initial_stride = i_stride;
        %                                 correction_x_toilet = positions(initial_stride, 1) - 9;
        %                                 correction_y_toilet = positions(initial_stride, 2) - 8.2;
        %                                 positions(initial_stride:end,1) = positions(initial_stride:end,1) -correction_x_toilet;
        %                                 positions(initial_stride:end,2) = positions(initial_stride:end,2) -correction_y_toilet;
        %                                 x_prev = 9;
        %                                 y_prev = 8.2;
        %                                 no_movement = false;
        %                                 toilet_already_detected = true;
        %                                 break
        %                             end
        %
        %                         case 3 % broom
        %                              if ~broom_already_detected
        %                                 initial_stride = i_stride;
        %                                 correction_x_broom = positions(initial_stride, 1) - 7.2;
        %                                 correction_y_broom = positions(initial_stride, 2) - 8.2;
        %                                 positions(initial_stride:end,1) = positions(initial_stride:end,1) -correction_x_broom;
        %                                 positions(initial_stride:end,2) = positions(initial_stride:end,2) -correction_y_broom;
        %                                 x_prev = 7.2;
        %                                 y_prev = 8.2;
        %                                 no_movement = false;
        %                                 broom_already_detected = true;
        %                                 break
        %                             end
        %
        %                         case 4 % pitcher
        %
        %                             if ~pitcher_already_detected
        %                                 initial_stride = i_stride;
        %                                 correction_x_pitcher = positions(initial_stride, 1) - 13.4;
        %                                 correction_y_pitcher = positions(initial_stride, 2) - 7.3;
        %                                 positions(initial_stride:end,1) = positions(initial_stride:end,1) - correction_x_pitcher;
        %                                 positions(initial_stride:end,2) = positions(initial_stride:end,2) - correction_y_pitcher;
        %                                 x_prev = 13.4;
        %                                 y_prev = 7.3;
        %                                 no_movement = false;
        %                                 pitcher_already_detected = true;
        %                                 break
        %                             end
        %
        %                         case 5 % brush
        %
        %                             if ~brush_already_detected
        %                                 initial_stride = i_stride;
        %                                 correction_x_brush = positions(initial_stride, 1) - 9;
        %                                 correction_y_brush = positions(initial_stride, 2) - 6.3;
        %                                 positions(initial_stride:end,1) = positions(initial_stride:end,1) -correction_x_brush;
        %                                 positions(initial_stride:end,2) = positions(initial_stride:end,2) -correction_y_brush;
        %                                 x_prev = 9;
        %                                 y_prev = 6.3;
        %                                 no_movement = false;
        %                                 brush_already_detected = true;
        %                                 break
        %                             end
        %                     end
        %                 end
        %             end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
end




% for i_lm=1:N_LM
%     hold on
%     plot(particles(1,1).Lm(i_lm,1),particles(1,1).Lm(i_lm,2),'d');
% end

function [particles, particles_weights, trajectories, current_best_particle] = ...
    update_resample(particles, rssi, j_rssi, particles_weights, N_PARTICLES, trajectories, THRESHOLD_RESAMPLE, previous_best_particle)
    
    % To update
    [particles, particles_weights, current_best_particle] = update(particles, rssi, j_rssi, particles_weights, N_PARTICLES, previous_best_particle);

    % To resample
    neff = 1 / sum(sqrt(particles_weights));
    if neff < THRESHOLD_RESAMPLE
        %indexes = double(py.mifunc.stratified_resample(particles_weights))+1;
        indexes = double(py.mifunc.systematic_resample(particles_weights))+1;
        [particles, particles_weights, trajectories] = resample_from_index(particles, particles_weights, indexes, N_PARTICLES, trajectories);
        %particles_weights = ones(1,N_PARTICLES) / N_PARTICLES; % PARTICLES WEIGHTS ARE INITILIZED TO AND EQUAL VALUE

    end
end

function [particles, particles_weights, current_best_particle] = update(particles, rssi, lm_id, particles_weights, N_PARTICLES, previous_best_particle)
    
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
    best_particles = [];
    biggest_weight = -1;
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).W = particles_weights(i_particle);
        
        % Get the best particle (particle with the biggest weight)
        if particles_weights(i_particle) >= biggest_weight
            best_particles = [best_particles, i_particle];
        end
        if find(best_particles==previous_best_particle)
            current_best_particle = previous_best_particle;
        else
            current_best_particle = best_particles(1); % for reproducibility the first item is always chosen.
        end
        
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

    end
    pause(0.05)
end











    