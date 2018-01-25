clear all;
close all;
clc
format long

delta_val = 1e-6;
stop = 200;

AAM_FILE = 'ESMGFZ_AAM_v1.0_03h_2004.asc';
HAM_FILE = 'ESMGFZ_HAM_v1.2_24h_2004.asc';
OAM_FILE = 'ESMGFZ_OAM_v1.0_03h_2004.asc';

FILE = OAM_FILE;
numParam    = 1:6; % Set from - to, to calculate for the parameters inside the inverall. See table below
maxNumParam = 6;

% Params: 
% 1: influence of HAM
% 2: influence of AAM
% 3: influence of OAM
% 4: trace
% 5: k_re
% 6: k_im

%% Creating the Text files used to read the w_0 and w_deltas

% delete all previous files
[timespan, step, h, grav_potent_ham, r_moon,r_sun, reference] = read_data(HAM_FILE);
[timespan, step, h, grav_potent_oam, r_moon,r_sun, reference] = read_data(AAM_FILE);
[timespan, step, h, grav_potent_aam, r_moon,r_sun, reference] = read_data(OAM_FILE);

% grav_potent_oam = zeros(5,96433);
% grav_potent_aam = zeros(5,96433);

initial = reference(:,1) .*3600;              % [rad / h]
G       = (6.674e-11) * (3600 * 3600);        % [m^3/(kg* h^2)]
GM_sun  = (1.32712442076e20) * (3600 * 3600); % [m^3/ h^2]
GM_moon = (4.9027779e12) * (3600 * 3600);     % [m^3/ h^2]
omega_N = (7.2921151467064e-5) * 3600;        % [rad/h]
Mass    = 5.9737e24;                          % [kg]
R       = 6378136.6;                          % [m]
A       = 0.3296108 * Mass * R * R;           % [kg * m^2]
B       = 0.3296108 * Mass * R * R;           % [kg * m^2]
C       = 0.3307007 * Mass * R * R;           % [kg * m^2]

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

h(:,end+1) = h(:,end);

disp('Delete all omega files before creating new ones')
for i = 1:maxNumParam
    if exist(['w_', num2str(i),'.txt'], 'file')==2
      disp(['Deleting: w_', num2str(i),'.txt'])
      delete(['w_', num2str(i),'.txt']);
    end
end

disp('Creating Files to read in omegas')
    
for i = 1:(max(numParam) - min(numParam) +2)
    disp(['creating file: w_', num2str(i) ,'.txt'])
    
    delta    = zeros(1,maxNumParam);
    delta(i) = delta_val;
    if i == (max(numParam) - min(numParam) +2)
        delta(i) = 0;
    end
        
    % All parameters influenced by delta
    c_20 = (1 + delta(1)) .* grav_potent_ham(1,:) + ...
           (1 + delta(2)) .* grav_potent_aam(1,:) + ...
           (1 + delta(3)) .* grav_potent_oam(1,:);
       
    c_21 = (1 + delta(1)) .* grav_potent_ham(2,:) + ...
           (1 + delta(2)) .* grav_potent_aam(2,:) + ...
           (1 + delta(3)) .* grav_potent_oam(2,:);
       
    s_21 = (1 + delta(1)) .* grav_potent_ham(3,:) + ...
           (1 + delta(2)) .* grav_potent_aam(3,:) + ...
           (1 + delta(3)) .* grav_potent_oam(3,:);
       
    c_22 = (1 + delta(1)) .* grav_potent_ham(4,:) + ...
           (1 + delta(2)) .* grav_potent_aam(4,:) + ...
           (1 + delta(3)) .* grav_potent_oam(4,:);
       
    s_22 = (1 + delta(1)) .* grav_potent_ham(5,:) + ...
           (1 + delta(2)) .* grav_potent_aam(5,:) + ...
           (1 + delta(3)) .* grav_potent_oam(5,:);
           
    tr      = (A + B + C) * (1 + delta(4));       % [kg * m^2]
    k_re    = 0.3077 * (1 + delta(5));            % [-] 
    k_im    = 0.0036 * (1 + delta(6));            % [-] 
    
    % calculating Tensors

    T_g_c1 = coefficient_T_g .* [((1/sqrt(3)) .* c_20 - c_22); -s_22; -c_21];
    T_g_c2 = coefficient_T_g .* [-s_22; ((1/sqrt(3)) .* c_20 + c_22); -s_21];
    T_g_c3 = coefficient_T_g .* [-c_21; -s_21; -(2/sqrt(3)) .* c_20];

    T_g_c1(1,:) = T_g_c1(1,:) + tr/3;
    T_g_c2(2,:) = T_g_c2(2,:) + tr/3;
    T_g_c3(3,:) = T_g_c3(3,:) + tr/3;

    T_g_c1(:,end+1) = T_g_c1(:,end);
    T_g_c2(:,end+1) = T_g_c2(:,end);
    T_g_c3(:,end+1) = T_g_c3(:,end);

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
    w = initial;
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

% for plotting
plot_w = [];

for i = 1:3:numel(w_delta(:,end))
    plot_w = [plot_w, [w_delta(i,end); w_delta(i+1,end); w_delta(i+2,end)]];
end

xp        = (R/omega_N) .* plot_w(1,:);
yp        = (R/omega_N) .* plot_w(2,:);

plot(xp, yp)
%% Calculation

disp('Creating A-Matrix')
diff = [];
for i = 1:(max(numParam) - min(numParam) +1)
    diff = [diff , [w_delta(:,i) - w_0]];
    A_mat(:,i) = (w_delta(:,i) - w_0) ./ delta_val;
end

disp('Calculating delta_x_dach')
% Singulärwertzerlegung
[U,S,V] = svd(A_mat,'econ');
    
% Singulärwertmatrix
lambdaInv = inv(S);
    
% Berechnung der Parameter herkoemmlich vs SVD
    
pseudo_inverse = V * lambdaInv * U';
% pseudo_inverse = ((A_mat' * A_mat) \ A_mat');
delta_x_dach =  pseudo_inverse * delta_l;



