clear all;
close all;
clc
format long

% Logging all console output
if exist('LogOutput.txt', 'file')==2
          disp('Deleting: LogOutput.txt')
          delete('LogOutput.txt');
end
diary('LogOutput.txt')
diary on

delta_val = 1e-6;
stop = 365;
stop_iterating = 1;
iteration = 1;

AAM_FILE = 'ESMGFZ_AAM_v1.0_03h_2004.asc';
HAM_FILE = 'ESMGFZ_HAM_v1.2_24h_2004.asc';
OAM_FILE = 'ESMGFZ_OAM_v1.0_03h_2004.asc';

numParam    = 1:6; % Set from - to, to calculate for the parameters inside the inverall. See table below
maxNumParam = 6;

% Params: 
% 1: influence of HAM
% 2: influence of AAM
% 3: influence of OAM
% 4: trace
% 5: k_re
% 6: k_im

x_vec = ones(maxNumParam,1) .* delta_val;

%% Read in all relevant data

[timespan, step, h, grav_potent_ham, r_moon,r_sun, reference] = read_data(HAM_FILE);
[timespan, step, h, grav_potent_oam, r_moon,r_sun, reference] = read_data(AAM_FILE);
[timespan, step, h, grav_potent_aam, r_moon,r_sun, reference] = read_data(OAM_FILE);

% grav_potent_oam = zeros(5,96433);
% grav_potent_aam = zeros(5,96433);

% Values that are not dependend on x, all depndend values are calculated in
% the calculate_w function
[initial, G, GM_sun, GM_moon, omega_N,...
 Mass, R, A, B, C, coefficient_F,...
 coefficient_T_g, coefficient_T_r, M] = undepentend_values(r_moon, r_sun, reference);

%% Main Loop
while iteration <= stop_iterating
    % delete all previous files
    disp('Delete all omega files before creating new ones')
    for i = 1:maxNumParam
        if exist(['w_', num2str(i),'.txt'], 'file')==2
          disp(['Deleting: w_', num2str(i),'.txt'])
          delete(['w_', num2str(i),'.txt']);
        end
    end

    disp('Creating Files to read in omegas')
    % creating the files to read in 
    for i = 1:(max(numParam) - min(numParam) +2)
        disp(['creating file: w_', num2str(i) ,'.txt'])

        x    = zeros(1,maxNumParam);
        if i == (max(numParam) - min(numParam) +2)
            x(i) = 0;
        else
            x(i) = x_vec(i);
        end

        result_m = calcuate_w(x, initial, M, h, grav_potent_ham,...
                              grav_potent_aam, grav_potent_oam, A, B, C, ...
                              coefficient_T_g, coefficient_T_r, coefficient_F, ...
                              timespan, stop);

        writetable(table(result_m), ['w_' , num2str(i), '.txt'])

        % Clearing aaray to avoid conflicts
        T_g_c1 = [];
        T_g_c2 = [];
        T_g_c3 = [];
    end

    %% Reading in w_0 and w_deltas
    disp('Reading in omegas and calculating delta_l')
    [timespan, step, h, grav_potent,r_moon,r_sun, reference] = read_data(HAM_FILE);

    w_0_struct = importdata(['w_', num2str(max(numParam) - min(numParam) +2), '.txt']);
    w_0        = w_0_struct.data;

    for i = 1:(max(numParam) - min(numParam) +1)
        disp(['Read in: w_', num2str(i), '.txt'])
        w_delta_struct = importdata(['w_', num2str(i), '.txt']);
        w_delta(:,i)   = w_delta_struct.data;
    end

    l_x = reference(1, 1:(numel(w_0)/3));
    l_y = reference(2, 1:(numel(w_0)/3));
    l_z = reference(3, 1:(numel(w_0)/3));
    l   = [];

    for i = 1:numel(l_x)
        l = [l, l_x(i), l_y(i), l_z(i)];
    end

    l = l';

    delta_l = l - w_0;

    %% Calculation of the A-matrix
    disp('Creating A-Matrix')
    diff = [];
    for i = 1:(max(numParam) - min(numParam) +1)
        diff = [diff , [w_delta(:,i) - w_0]];
        A_mat(:,i) = (w_delta(:,i) - w_0) ./ delta_val;
    end
    
    %% Calculatin the new x-vec
    disp('Calculating delta_x_dach')
    % Singulärwertzerlegung
    [U,S,V] = svd(A_mat,'econ');

    % Singulärwertmatrix
    lambdaInv = inv(S);

    % Berechnung der Parameter herkoemmlich vs SVD

    pseudo_inverse = V * lambdaInv * U';
    % pseudo_inverse = ((A_mat' * A_mat) \ A_mat');
    delta_x_dach =  pseudo_inverse * delta_l

    x_dach = x_vec + delta_x_dach

    x_vec = x_dach;

    iteration = iteration +1;

end 

%% Plot the final results
result_m = calcuate_w(x_dach, initial, M, h, grav_potent_ham,...
                      grav_potent_aam, grav_potent_oam, A, B, C, ...
                      coefficient_T_g, coefficient_T_r, coefficient_F, ...
                      timespan, stop);
                  
                 

% prepare the values for plotting
plot_w = [];

for i = 1:3:numel(result_m)
    plot_w = [plot_w, [result_m(i); result_m(i+1); result_m(i+2)]];
end

xp        = (R/omega_N) .* plot_w(1,:);
yp        = (R/omega_N) .* plot_w(2,:);

plot(xp, yp)

diary off