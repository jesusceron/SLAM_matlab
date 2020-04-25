clc;
clearvars -except positions step_events beac_rssi_fixed beac_rssi_fixed_filtered beac_rssi_activity beac_rssi_activity_filtered
close all;

%[positions, step_events, beac_rssi_fixed, beac_rssi_fixed_filtered, beac_rssi_activity, beac_rssi_activity_filtered] = Load_data();
%step_events = [1, step_events, length(beac_rssi_fixed)]; % Add '0' and the lenght of the dataset as points of support

rssi_room_unfiltered = beac_rssi_fixed(:,1);
rssi_kitchen_unfiltered = beac_rssi_fixed(:,2);
rssi_bathroom_unfiltered = beac_rssi_fixed(:,3);
rssi_dining_unfiltered = beac_rssi_fixed(:,4);
rssi_living_unfiltered = beac_rssi_fixed(:,5);
rssi_room = beac_rssi_fixed_filtered(:,1);
rssi_kitchen = beac_rssi_fixed_filtered(:,2);
rssi_bathroom = beac_rssi_fixed_filtered(:,3);
rssi_dining = beac_rssi_fixed_filtered(:,4);
rssi_living = beac_rssi_fixed_filtered(:,5);

rssi_door_unfiltered = beac_rssi_activity(:,1);
rssi_toilet_unfiltered = beac_rssi_activity(:,2);
rssi_broom_unfiltered = beac_rssi_activity(:,3);
rssi_pitcher_unfiltered = beac_rssi_activity(:,4);
rssi_brush_unfiltered = beac_rssi_activity(:,5);
rssi_door = beac_rssi_activity_filtered(:,1);
rssi_toilet = beac_rssi_activity_filtered(:,2);
rssi_broom = beac_rssi_activity_filtered(:,3);
rssi_pitcher = beac_rssi_activity_filtered(:,4);
rssi_brush = beac_rssi_activity_filtered(:,5);

beacon_fixed = [1,1,1,1,1,1,1,1,1,1];

%% CREATION OF PARTICLES
N_PARTICLES = 100;
N_LM = 8;
LM_SIZE = 2;
THRESHOLD_RESAMPLE = N_PARTICLES / 10;
particles = [];
particles_weights = ones(1,N_PARTICLES) / N_PARTICLES;

STD_PERSON_POSITION = 0.05;

for i_particle=1:N_PARTICLES
    particles = [particles, Particle(N_PARTICLES,N_LM,LM_SIZE)];
end
particles = init_beacons_position(particles,N_PARTICLES);

%% MAIN LOOP
x_prev = 0;
y_prev = 0;
for i_stride=1:length(positions)
    display(i_stride)
    % PREDICTION
    x = positions(i_stride, 1);
    y = positions(i_stride, 2);
    for i_particle=1:N_PARTICLES
        dist_x = (x - x_prev) + randn() * STD_PERSON_POSITION;
        dist_y = (y - y_prev) + randn() * STD_PERSON_POSITION;
        particles(1,i_particle).X = particles(i_particle).X + dist_x;
        particles(1,i_particle).Y = particles(i_particle).Y + dist_y;
        
        particles(1,i_particle).T_x = [particles(1,i_particle).T_x, particles(1,i_particle).X];
        particles(1,i_particle).T_y = [particles(1,i_particle).T_y, particles(1,i_particle).Y];
    end    
    x_prev = x;
    y_prev = y;
    
    initial_sample = step_events(i_stride);
    final_sample = step_events(i_stride + 1) - 1;
    % UPDATE
    rssis = [rssi_room(initial_sample:final_sample), rssi_kitchen(initial_sample:final_sample),...
            rssi_bathroom(initial_sample:final_sample), rssi_dining(initial_sample:final_sample),...
            rssi_living(initial_sample:final_sample), rssi_door(initial_sample:final_sample),...
            rssi_toilet(initial_sample:final_sample), rssi_brush(initial_sample:final_sample)];

    for i_rssis=1:length(rssis)
        non_zero_indexes = find(rssis(i_rssis));
        if ~isempty(non_zero_indexes)
        
            for j_rssi=non_zero_indexes
                rssi = rssis(i_rssis,j_rssi);
                if (rssi ~=0) && (rssi > -88)

                    [particles, particles_weights] = update(particles, rssi, j_rssi, particles_weights, N_PARTICLES);

                    % RESAMPLE
                    neff = 1 / sum(sqrt(particles_weights));
                    if neff < THRESHOLD_RESAMPLE
                        indexes = double(py.mifunc.stratified_resample(particles_weights))+1;
                        %indexes = py.mifunc.systematic_resample(particles_weights);
                        particles = resample_from_index(particles, indexes, N_PARTICLES);

                        particles_weights = ones(1,N_PARTICLES) / N_PARTICLES; % PARTICLES WEIGHTS ARE INITILIZED TO AND EQUAL VALUE
                    end

                end
            end
        end
    end
        
end

for i_particle=1:N_PARTICLES
    plot(particles(1,i_particle).T_x, particles(1,i_particle).T_y);
    hold on
end

for i_lm=1:N_LM
    plot(particles(1,1).Lm(i_lm,1),particles(1,1).Lm(i_lm,2),'d');
end


function [particles, particles_weights] = update(particles, rssi, lm_id, particles_weights, N_PARTICLES)
    
    distance_rssi = 10 ^ ((rssi + 83) / -10 * 0.6);  % Calculate distance to beacon
    R = distance_rssi / 10;

    for i_particle=1:N_PARTICLES
        % Get the distance between the person and the landmark
        dx = particles(1,i_particle).X - particles(1,i_particle).Lm(lm_id,1);
        dy = particles(i_particle).Y - particles(i_particle).Lm(lm_id,2);
        dist = (dx * dx) + (dy * dy);
        residual = distance_rssi - dist;

        % Compute Jacobians
        H = [-dx / dist, -dy / dist];

        % Compute covariance of the residual
        % covV = H * Cov_s * H^T + error
        HxCov = [particles(1,i_particle).LmP{1,lm_id}(1, 1) * H(1) + particles(1,i_particle).LmP{1,lm_id}(1, 2) * H(2),...
                 particles(1,i_particle).LmP{1,lm_id}(2, 1) * H(1) + particles(1,i_particle).LmP{1,lm_id}(2, 2) * H(2)];

        covV = (HxCov(1) * H(1)) + (HxCov(2) * H(2)) + R;

        % Calculate Kalman gain
        K_gain = [HxCov(1) * (1 / covV), HxCov(1) * (1.0 / covV)];

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
    end
    particles_weights = particles_weights / sum(particles_weights);
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).W = particles_weights(i_particle);
    end
end

function particles = resample_from_index(particles, indexes, N_PARTICLES)
    % particles_weights.fill(1.0 / N_PARTICLES)
    for i_particle = fliplr(linspace(1,N_PARTICLES,N_PARTICLES))
        particles(1,i_particle).W = 1 / N_PARTICLES;
        particles(1,i_particle).X = particles(1,indexes(i_particle)).X;
        particles(1,i_particle).Y = particles(1,indexes(i_particle)).Y;
        particles(1,i_particle).Lm = particles(1,indexes(i_particle)).Lm;
        particles(1,i_particle).LmP = particles(1,indexes(i_particle)).LmP;
        particles(1,i_particle).T_x = particles(1,indexes(i_particle)).T_x;
        particles(1,i_particle).T_y = particles(1,indexes(i_particle)).T_y;
    end
end











    