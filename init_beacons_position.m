function particles = init_beacons_position(particles,N_PARTICLES)
% Initialization of beacons position

    % Initialize Lm[0] (beacon in the room)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(1,:) = [16.5, 6.3];
        particles(1,i_particle).LmP{1} = [[0.00287, 0]; [0, 0.00061]];
    end

    % Initialize Lm[1] (beacon in the kitchen)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(2,:) = [11.6, 7.9];
        particles(1,i_particle).LmP{2} = [[0.00287, 0]; [0, 0.00061]];
    end

    % Initialize Lm[2] (beacon in the bathroom)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(3,:) = [8, 8.5];
        particles(1,i_particle).LmP{3} = [[0.00287, 0]; [0, 0.00061]];
    end

    % Initialize Lm[3] (beacon in the dining room)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(4,:) = [12.1, 2.8];
        particles(1,i_particle).LmP{4} = [[0.00287, 0]; [0, 0.00061]];
    end

    % Initialize Lm[4] (beacon in the living room)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(5,:) = [15.7, 2.4];
        particles(1,i_particle).LmP{5} = [[0.00287, 0]; [0, 0.00061]];
    end

    % Initialize Lm[4] (beacon in the door)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(6,:) = [5.4, 5.5];
        particles(1,i_particle).LmP{6}= [[0.00287, 0]; [0, 0.00061]];
    end

    % Initialize Lm[4] (beacon in the toilet)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(7,:) = [9, 8.2];
        particles(1,i_particle).LmP{7} = [[0.00287, 0]; [0, 0.00061]];
    end

    % % Initialize Lm[4] (beacon in the broom)
    % for i_particle=1:N_PARTICLES
    %     particles(1,i_particle).Lm[4][:] = [7.2, 8.2]
    %     particles(1,i_particle).LmP[4][:] = [[0.00287, 0], [0, 0.00061]]

    % % Initialize Lm[4] (beacon in the pitcher)
    % for i_particle=1:N_PARTICLES
    %     particles(1,i_particle).Lm[4][:] = [13.4, 7.3]
    %     particles(1,i_particle).LmP[4][:] = [[0.00287, 0]; [0, 0.00061]]

    % Initialize Lm[4] (beacon in the brush)
    for i_particle=1:N_PARTICLES
        particles(1,i_particle).Lm(8,:) = [9, 6.3];
        particles(1,i_particle).LmP{8} = [[0.00287, 0]; [0, 0.00061]];
    end
end