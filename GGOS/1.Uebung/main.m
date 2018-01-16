clear all;
close all;
clc
%% Feeding in data
AAM_FILE = 'ESMGFZ_AAM_v1.0_03h_2004.asc';
HAM_FILE = 'ESMGFZ_HAM_v1.2_24h_2004.asc';
OAM_FILE = 'ESMGFZ_OAM_v1.0_03h_2004.asc';

% To get another set of data for the calculation change the file names. The
% sample rate will be adapted automatically
[timespan, step, h, grav_potent,r_moon,r_sun, reference] = read_data(OAM_FILE);
r_moon = r_moon;
r_sun = r_sun;
initial = reference(:,1);
c_20 = grav_potent(1, :);
c_21 = grav_potent(2, :);
s_21 = grav_potent(3, :);
c_22 = grav_potent(4, :);
s_22 = grav_potent(5, :);

%% General Values

G       = (6.674e-11) * (3600 * 3600);        % [m^3/(kg* h^2)]
GM_sun  = (1.32712442076e20) * (3600 * 3600); % [m^3/ h^2]
GM_moon = (4.9027779e12) * (3600 * 3600);     % [m^3/ h^2]
omega_N = (7.2921151467064e-5) * 3600;        % [rad/h]
Mass    = 5.9737e24;                          % [kg]
R       = 6378136.6;                          % [m]
A       = 0.3296108 * Mass * R * R;           % [kg * m^2]
B       = 0.3296108 * Mass * R * R;           % [kg * m^2]
C       = 0.3307007 * Mass * R * R;           % [kg * m^2]
tr      = A + B + C;                          % [kg * m^2]
k_re    = 0.3077;                             % [-] 
k_im    = 0.0036;                             % [-] 

coefficient_T_g = sqrt(5/3) * Mass * R * R;
coefficient_T_r = (omega_N * R^5) / (3 * G);
coefficient_F   = (omega_N * omega_N * R^5) / (3 * G);

%%
dist_sun   = sqrt(r_sun(1,:).^2 + r_sun(2,:).^2 + r_sun(3,:).^2);
dist_moon   = sqrt(r_moon(1,:).^2 + r_moon(2,:).^2 + r_moon(3,:).^2);

M_sun(1,:) = (((C-B) .* r_sun(2,:) .* r_sun(3,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);
M_sun(2,:) = (((A-C) .* r_sun(1,:) .* r_sun(3,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);
M_sun(3,:) = (((B-A) .* r_sun(1,:) .* r_sun(2,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);

M_moon(1,:) = (((C-B) .* r_moon(2,:) .* r_moon(3,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
M_moon(2,:) = (((A-C) .* r_moon(1,:) .* r_moon(3,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
M_moon(3,:) = (((B-A) .* r_moon(1,:) .* r_moon(2,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);

M = M_sun + M_moon;

%% 
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


%% Calculation
dt = 1;
counter = 1;
w = initial.*3600;
for t = 1:dt:24*365
    
    w_x = w(1);
    w_y = w(2);
    w_z = w(3);
    
    T_r = coefficient_T_r .* [0,     0,      k_re * w_x + k_im * w_y;
                              0,     0,      k_re * w_y - k_im * w_x;
            k_re * w_x + k_im * w_y,     k_re * w_y - k_im * w_x,         0];
        
    disp ('Tag: ')
    disp (t/24) 
    F = T_g{t} + T_r + coefficient_F .* [k_re, k_im, 0;
                                        -k_im, k_re, 0;
                                         0   ,  0  , 0];
    
    w = w + dt .* (F\(M(:,t) - (dT_g{t} * w) - cross(w, (T_g{t} + T_r) * w) - cross(w,h(:,t)) - dh{t})); 
    
    result{t} = w;
end
%%
result_m = cell2mat(result);

xp        = (R/omega_N) .* result_m(1,:);
yp        = (R/omega_N) .* result_m(2,:);

xp_reference = (R/omega_N) .* reference(1,:) .* 3600;
yp_reference = (R/omega_N) .* reference(2,:) .* 3600;

delta_LOD = 86400/omega_N .* (omega_N - result_m(3,:)./3600);
LOD = 8.6376e+04 - delta_LOD;
LOD = LOD .* 1000;
figure(1)
plot(xp,yp)
title('Polar Motion at Earth Surface')
ylabel('y[m]')
xlabel('x[m]')

figure(2)
plot(delta_LOD)
title('Length of day')
ylabel('Length of Day [s]')
xlabel('Time [h]')

figure(3)
subplot(3,1,1)
plot(result_m(1,:)./3600)
title('Angular velocity - w_x')
ylabel('w_x [rad/s]')
xlabel('Time [h]')
subplot(3,1,2)
plot(result_m(2,:)./3600)
title('Angular velocity - w_y')
ylabel('w_y [rad/s]')
xlabel('Time [h]')
subplot(3,1,3)
plot(result_m(3,:)./3600)
title('Angular velocity - w_z')
ylabel('w_z [rad/s]')
xlabel('Time [h]')

figure(4)
plot(LOD)
title('Delta Length of day')
ylabel('Length of Day [ms]')
xlabel('Time [h]')






