function particles = init_beacons_position(particles,N_PARTICLES)
% Initialization of beacons position

%     % Initialize Lm[5] (beacon in the door)
%     for i_particle=1:N_PARTICLES
%         particles(1,i_particle).Lm(1,:) = [5, 5.5];
%         particles(1,i_particle).LmP{1}= [[0, 0]; [0, 0]];
%     end

%     % Initialize Lm[6] (beacon in the toilet)
%     for i_particle=1:N_PARTICLES
%         particles(1,i_particle).Lm(2,:) = [9.1, 8.3];
%         particles(1,i_particle).LmP{2} = [[0, 0]; [0, 0]];
%     end

%     % Initialize Lm[7] (beacon in the broom)
%     for i_particle=1:N_PARTICLES
%         particles(1,i_particle).Lm(3,:) = [7.3, 8];
%         particles(1,i_particle).LmP{3} = [[0, 0]; [0, 0]];
%     end

%     % Initialize Lm[8] (beacon in the pitcher)
%     for i_particle=1:N_PARTICLES
%         particles(1,i_particle).Lm(4,:) = [13.5, 7.4];
%         particles(1,i_particle).LmP{4} = [[0, 0]; [0, 0]];
%     end

%     % Initialize Lm[9] (beacon in the brush)
%     for i_particle=1:N_PARTICLES
%         particles(1,i_particle).Lm(5,:) = [8.7, 6.3];
%         particles(1,i_particle).LmP{5} = [[0, 0]; [0, 0]];
%     end
    % Initialize Lm[0] (beacon in the room)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(1,:) = [16.6, 6.4];
        particles(1,i_particle).LmP{1} = [[0, 0]; [0, 0]];
    end

    % Initialize Lm[1] (beacon in the kitchen)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(2,:) = [11.7, 8];
        particles(1,i_particle).LmP{2} = [[0, 0]; [0, 0]];
    end

    % Initialize Lm[2] (beacon in the bathroom)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(3,:) = [8.1, 8.6];
        particles(1,i_particle).LmP{3} = [[0, 0]; [0, 0]];
    end

    % Initialize Lm[3] (beacon in the dining room)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(4,:) = [12.2, 3];
        particles(1,i_particle).LmP{4} = [[0, 0]; [0, 0]];
    end

    % Initialize Lm[4] (beacon in the living room)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(5,:) = [15.8, 2.5];
        particles(1,i_particle).LmP{5} = [[0, 0]; [0, 0]];
    end


end