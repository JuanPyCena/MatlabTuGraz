clear all;
close all;
clc
format long

delta_val = 1e-15;
stop = 365;

AAM_FILE = 'ESMGFZ_AAM_v1.0_03h_2004.asc';
HAM_FILE = 'ESMGFZ_HAM_v1.2_24h_2004.asc';
OAM_FILE = 'ESMGFZ_OAM_v1.0_03h_2004.asc';

FILE = OAM_FILE;
numParam    = 2;
maxNumParam = 10;
%% Creating the Text files used to read the w_0 and w_deltas

disp('Creating Files to read in omegas')
[timespan, step, h, grav_potent,r_moon,r_sun, reference] = read_data(FILE);

for i = 1:(numParam+1)
    disp(['w_', num2str(i) ,'.txt'])
    delta    = zeros(1,maxNumParam);
    delta(i) = delta_val;
    if i == numParam+1
        delta(i) = 0;
    end
    
    r_moon = r_moon;
    r_sun = r_sun;
    initial = reference(:,1);
    c_20 = grav_potent(1, :) + delta(3);
    c_21 = grav_potent(2, :) + delta(4);
    s_21 = grav_potent(3, :) + delta(5);
    c_22 = grav_potent(4, :) + delta(6);
    s_22 = grav_potent(5, :) + delta(7);

    G       = (6.674e-11) * (3600 * 3600);        % [m^3/(kg* h^2)]
    GM_sun  = (1.32712442076e20) * (3600 * 3600); % [m^3/ h^2]
    GM_moon = (4.9027779e12) * (3600 * 3600);     % [m^3/ h^2]
    omega_N = (7.2921151467064e-5) * 3600;        % [rad/h]
    Mass    = 5.9737e24;                          % [kg]
    R       = 6378136.6;                          % [m]
    A       = 0.3296108 * Mass * R * R + delta(8);           % [kg * m^2]
    B       = 0.3296108 * Mass * R * R + delta(9);           % [kg * m^2]
    C       = 0.3307007 * Mass * R * R + delta(10);           % [kg * m^2]
    tr      = A + B + C;                          % [kg * m^2]
    k_re    = 0.3077 + delta(1);                             % [-] 
    k_im    = 0.0036 + delta(2);                             % [-] 

    coefficient_T_g = sqrt(5/3) * Mass * R * R;
    coefficient_T_r = (omega_N * R^5) / (3 * G);
    coefficient_F   = (omega_N * omega_N * R^5) / (3 * G);

    dist_sun   = sqrt(r_sun(1,:).^2 + r_sun(2,:).^2 + r_sun(3,:).^2);
    dist_moon   = sqrt(r_moon(1,:).^2 + r_moon(2,:).^2 + r_moon(3,:).^2);

    M_sun(1,:) = (((C-B) .* r_sun(2,:) .* r_sun(3,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);
    M_sun(2,:) = (((A-C) .* r_sun(1,:) .* r_sun(3,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);
    M_sun(3,:) = (((B-A) .* r_sun(1,:) .* r_sun(2,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);

    M_moon(1,:) = (((C-B) .* r_moon(2,:) .* r_moon(3,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
    M_moon(2,:) = (((A-C) .* r_moon(1,:) .* r_moon(3,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
    M_moon(3,:) = (((B-A) .* r_moon(1,:) .* r_moon(2,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);

    M = M_sun + M_moon;
    T_g_c1 = coefficient_T_g .* [((1/sqrt(3)) .* c_20 - c_22); -s_22; -c_21];
    T_g_c2 = coefficient_T_g .* [-s_22; ((1/sqrt(3)) .* c_20 + c_22); -s_21];
    T_g_c3 = coefficient_T_g .* [-c_21; -s_21; -(2/sqrt(3)) .* c_20];

    T_g_c1(1,:) = T_g_c1(1,:) + tr/3;
    T_g_c2(2,:) = T_g_c2(2,:) + tr/3;
    T_g_c3(3,:) = T_g_c3(3,:) + tr/3;

    T_g_c1(:,end+1) = T_g_c1(:,end);
    T_g_c2(:,end+1) = T_g_c2(:,end);
    T_g_c3(:,end+1) = T_g_c3(:,end);


    h(:,end+1) = h(:,end);

    for time = 1:timespan+1
        T_g{time}  = [T_g_c1(:, time),T_g_c2(:, time), T_g_c3(:, time)];
    end

    for time = 1:timespan
        dT_g{time} = T_g{time + 1} - T_g{time};
    end

    for time = 1:(length(h)-1)
        dh{time} = h(:,time + 1) - h(:,time);
    end

    dt = 1;
    counter = 1;
    w = initial.*3600;
    result_m = [];
    for t = 1:dt:24*stop

        w_x = w(1);
        w_y = w(2);
        w_z = w(3);

        T_r = coefficient_T_r .* [0,     0,      k_re * w_x + k_im * w_y;
                                  0,     0,      k_re * w_y - k_im * w_x;
                k_re * w_x + k_im * w_y,     k_re * w_y - k_im * w_x,         0];

        F = T_g{t} + T_r + coefficient_F .* [k_re, k_im, 0;
                                            -k_im, k_re, 0;
                                             0   ,  0  , 0];

        w = w + dt .* (F\(M(:,t) - (dT_g{t} * w) - cross(w, (T_g{t} + T_r) * w) - cross(w,h(:,t)) - dh{t})); 

        result_m = [result_m, w(1), w(2), w(3)];
    end

    result_m = result_m';
    writetable(table(result_m), ['w_' , num2str(i), '.txt'])

end

%% Reading in w_0 and w_deltas
disp('Reading in omegas and calculating delta_l')
[timespan, step, h, grav_potent,r_moon,r_sun, reference] = read_data(FILE);

w_0_struct = importdata(['w_', num2str(numParam+1), '.txt']);
w_0        = w_0_struct.data;

for i = 1:numParam
    w_delta_struct = importdata(['w_', num2str(i), '.txt']);
    w_delta(:,i)   = w_delta_struct.data;
end

% w_delta_struct_1  = importdata('w_1.txt');
% w_delta_struct_2  = importdata('w_2.txt');
% w_delta_struct_3  = importdata('w_3.txt');
% w_delta_struct_4  = importdata('w_4.txt');
% w_delta_struct_5  = importdata('w_5.txt');
% w_delta_struct_6  = importdata('w_6.txt');
% w_delta_struct_7  = importdata('w_7.txt');
% w_delta_struct_8  = importdata('w_8.txt');
% w_delta_struct_9  = importdata('w_9.txt');
% w_delta_struct_10 = importdata('w_10.txt');
% 
% w_delta(:,1)  = w_delta_struct_1.data;
% w_delta(:,2)  = w_delta_struct_2.data;
% w_delta(:,3)  = w_delta_struct_3.data;
% w_delta(:,4)  = w_delta_struct_4.data;
% w_delta(:,5)  = w_delta_struct_5.data;
% w_delta(:,6)  = w_delta_struct_6.data;
% w_delta(:,7)  = w_delta_struct_7.data;
% w_delta(:,8)  = w_delta_struct_8.data;
% w_delta(:,9)  = w_delta_struct_9.data;
% w_delta(:,10) = w_delta_struct_10.data;

l_x = reference(1, 1:(numel(w_0)/3));
l_y = reference(2, 1:(numel(w_0)/3));
l_z = reference(3, 1:(numel(w_0)/3));
l   = [];

for i = 1:numel(l_x)
    l = [l, l_x(i), l_y(i), l_z(i)];
end

l = l';

delta_l    = l - w_0;

%% Calculation

disp('Creating A-Matrix')
for i = 1:numParam
    A_mat(:,i) = (w_delta(:,i) - w_0) ./ delta_val;
end

disp('Calculating delta_x_dach')
pseudo_inverse = ((A_mat' * A_mat) \ A_mat');
delta_x_dach =  pseudo_inverse * delta_l;



