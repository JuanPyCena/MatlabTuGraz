clear all;
close all;
clc
warning('off','all')
warning
format long

% Logging all console output
if exist('LogOutput.txt', 'file')==2
          disp('Deleting: LogOutput.txt')
          delete('LogOutput.txt');
end
diary('LogOutput.txt')
diary on

delta_val = 0
stop = 365
stop_iterating = 5
iteration = 1;
threshold = 100000

AAM_FILE = 'ESMGFZ_AAM_v1.0_03h_2004.asc';
HAM_FILE = 'ESMGFZ_HAM_v1.2_24h_2004.asc';
OAM_FILE = 'ESMGFZ_OAM_v1.0_03h_2004.asc';

numParam    = 1:6 % Set from - to, to calculate for the parameters inside the inverall. See table below
maxNumParam = 6

% Params: 
% 1: influence of HAM
% 2: influence of AAM
% 3: influence of OAM
% 4: trace
% 5: k_re
% 6: k_im

x_vec = ones(numel(numParam),1) .* delta_val

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
    disp(['Iteration: ', num2str(iteration)])
    disp('Delete all omega files before creating new ones')
    for i = 1:(maxNumParam+1)
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

        result_calculate = calcuate_w(x, initial, M, h, grav_potent_ham,...
                              grav_potent_aam, grav_potent_oam, A, B, C, ...
                              coefficient_T_g, coefficient_T_r, coefficient_F, ...
                              timespan, stop);
        T = table(result_calculate);
        writetable(T, ['w_' , num2str(i), '.txt'])

        % Clearing aaray to avoid conflicts
        T_g_c1 = [];
        T_g_c2 = [];
        T_g_c3 = [];
    end

    %% Reading in w_0 and w_deltas
    disp('Reading in omegas and calculating delta_l')
    [timespan, step, h, grav_potent,r_moon,r_sun, reference] = read_data(HAM_FILE);
    
    w_0 = table2array(readtable(['w_', num2str(max(numParam) - min(numParam) +2), '.txt']));

    for i = 1:(max(numParam) - min(numParam) +1)
        disp(['Read in: w_', num2str(i), '.txt'])
        w_delta(:,i) = table2array(readtable(['w_', num2str(i), '.txt']));
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
%         diff = [diff , [w_delta(:,i) - w_0]];
        A_mat(:,i) = (w_delta(:,i) - w_0) ./ delta_val;
    end
    
    %% Calculatin the new x-vec
    disp('Calculating delta_x_dach')
    
    % Preventing NAN and INF
    for i = 1:numel(A_mat)
       if isnan(A_mat(i)) | isinf(A_mat(i))
           A_mat(i) = 0;
       end
    end
    
    % Singulärwertzerlegung
    [U,S,V] = svd(A_mat,'econ');

    % Singulärwertmatrix
    lambdaInv = inv(S);

    % Berechnung der Parameter herkoemmlich vs SVD

    pseudo_inverse = V * lambdaInv * U';
    % pseudo_inverse = ((A_mat' * A_mat) \ A_mat');
    delta_x_dach   =  pseudo_inverse * delta_l

    x_dach    = x_vec + delta_x_dach
    x_vec     = x_dach;
    
    x_check = zeros(1,maxNumParam);

    for i = 1:numel(x_dach)
       x_check(i) = x_dach(i);
    end

    result_m = calcuate_w(x_check, initial, M, h, grav_potent_ham,...
                          grav_potent_aam, grav_potent_oam, A, B, C, ...
                          coefficient_T_g, coefficient_T_r, coefficient_F, ...
                          timespan, stop);
    
    % prepare the values for checking
    w_check = [];

    for i = 1:3:numel(result_m)
        w_check = [w_check, [result_m(i); result_m(i+1); result_m(i+2)]];
    end
    
    if (abs(max(w_check - reference(:, 1:length(w_check)))) < threshold)
        disp(['JOSEF: ', num2str(iteration)])
        break 
    end
    
    iteration = iteration +1;
    
    % Clear after every iteratio to avoid errors
    clear w_delta
    clear w_0
end 

%% Plot the final results
xp        = (R/omega_N) .* w_check(1,:);
yp        = (R/omega_N) .* w_check(2,:);

plot(xp, yp)

diary off