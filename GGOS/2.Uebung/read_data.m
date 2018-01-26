function [timespan, step, h, grav_potent,r_moon,r_sun, reference] = read_data(fileName_beginning, year_start, year_end)
% This function reads in a series of given files and then returns the
% timespan, motion terms of the angular momentum, gravitational potential,
% and the positions of the moon and sun. 

    % General constants
    omega_N = (7.2921151467064e-5) * 3600;        % [rad/s]
    Mass    = 5.9737e24;                          % [kg]
    R       = 6378136.6;                          % [m]
    A       = 0.3296108 * Mass * R * R;           % [kg * m^2]
    C       = 0.3307007 * Mass * R * R;           % [kg * m^2]

    % importing the data
    delimiterIn    = ' ';
    file_format    = '.asc';
    file_numbering = '';
    file_ending    = '';
    
    DATA_motionMomentum = [];
    for i = year_start:year_end
        if i < 10
            file_numbering = strcat('0', num2str(i));
        else
            file_numbering = num2str(i);
        end
        file_ending             = strcat(file_numbering, file_format);
        fileName_motionMomentum = strcat(fileName_beginning, file_ending);
        motion_data             = importdata(fileName_motionMomentum, delimiterIn);
        DATA_motionMomentum     = [DATA_motionMomentum; motion_data(:,9:11)]; 
    end
    
    DATA_gravPotential  = importdata('potentialCoefficientsAOHIS.txt', delimiterIn);
    DATA_moon           = importdata('moon.txt', delimiterIn);
    DATA_sun            = importdata('sun.txt', delimiterIn);
    DATA_reference      = importdata('earthRotationVector.txt', delimiterIn);
    
    % allocating the timespan
    timespan = length(DATA_gravPotential(:,1));
    
    % defining the timestep by checking the fileName
    if ~isempty(strfind(fileName_motionMomentum, '_24h_'))
        step = 24;
    elseif ~isempty(strfind(fileName_motionMomentum, '_03h_'))
        step = 3;
    else
        step = 1;
    end
    
    % allocating the motionMomentum
    h_read(1,:) = DATA_motionMomentum(:,1) .* (omega_N * (C - A) / (1.610));
    h_read(2,:) = DATA_motionMomentum(:,2) .* (omega_N * (C - A) / (1.610));
    h_read(3,:) = DATA_motionMomentum(:,3) .* (omega_N *  C / (1.125));
    
    h(1,:) = (interp1(h_read(1,:),1:1/step:length(h_read)))';
    h(2,:) = (interp1(h_read(2,:),1:1/step:length(h_read)))';
    h(3,:) = (interp1(h_read(3,:),1:1/step:length(h_read)))';

    % allocating the gravitational potential
    grav_potent(1,:) = DATA_gravPotential(1:timespan,2);
    grav_potent(2,:) = DATA_gravPotential(1:timespan,3);
    grav_potent(3,:) = DATA_gravPotential(1:timespan,4);
    grav_potent(4,:) = DATA_gravPotential(1:timespan,5);
    grav_potent(5,:) = DATA_gravPotential(1:timespan,6);
    
    % allocation of the moon's and the sun's coordinates
    r_moon(1,:) = DATA_moon(1:length(DATA_moon),2);
    r_moon(2,:) = DATA_moon(1:length(DATA_moon),3);
    r_moon(3,:) = DATA_moon(1:length(DATA_moon),4);
    
    r_sun(1,:) = DATA_sun(1:length(DATA_sun),2);
    r_sun(2,:) = DATA_sun(1:length(DATA_sun),3);
    r_sun(3,:) = DATA_sun(1:length(DATA_sun),4);
    
    % allocating the reference data
    reference(1,:) = DATA_reference(1:timespan,2);
    reference(2,:) = DATA_reference(1:timespan,3);
    reference(3,:) = DATA_reference(1:timespan,4);
end

