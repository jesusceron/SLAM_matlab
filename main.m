clc;
clearvars -except positions beac_rssi_fixed beac_rssi_fixed_filtered beac_rssi_activity beac_rssi_activity_filtered
close all;

%[positions, beac_rssi_fixed, beac_rssi_fixed_filtered, beac_rssi_activity, beac_rssi_activity_filtered] = Load_data();

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
N_PARTICLES = 1000;
N_LM = 5;
LM_SIZE = 2;
particles = [];

STD_PERSON_POSITION = 0.05;

for i_particle=1:N_PARTICLES
    particles = [particles, Particle(N_PARTICLES,N_LM,LM_SIZE)];
end

%% MAIN LOOP
x_prev = 0;
y_prev = 0;
for i_stride=1:length(positions)

    x = positions(i_stride, 1);
    y = positions(i_stride, 2);
    
    % PREDICTION
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
        
end

for i_particle=1:N_PARTICLES
    plot(particles(1,i_particle).T_x, particles(1,i_particle).T_y);
    hold on
end




    