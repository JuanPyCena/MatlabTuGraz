clear all;
close all;
clc
%% General Values
format long

AAM_FILE = 'ESMGFZ_AAM_v1.0_03h_2004.asc';
HAM_FILE = 'ESMGFZ_HAM_v1.2_24h_2004.asc';
OAM_FILE = 'ESMGFZ_OAM_v1.0_03h_2004.asc';

% To get another set of data for the calculation change the file names. The
% sample rate will be adapted automatically
[timespan, step, h, grav_potent,r_moon,r_sun, reference] = read_data(HAM_FILE);
w_initial = reference(:,1);
c_20_initial = grav_potent(1, 1);
c_21_initial = grav_potent(2, 1);
s_21_initial = grav_potent(3, 1);
c_22_initial = grav_potent(4, 1);
s_22_initial = grav_potent(5, 1);

G       = (6.674e-11);% * (3600 * 3600);        % [m^3/(kg* s^2)]
GM_sun  = (1.32712442076e20);% * (3600 * 3600); % [m^3/ s^2]
GM_moon = (4.9027779e12);% * (3600 * 3600);     % [m^3/ s^2]
omega_N = (7.2921151467064e-5);% * 3600;        % [rad/s]
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


%% Calculation
dt = 1;
w_in = w_initial;%.*3600                      % [m^3/ s^2]


x          = [c_20_initial, c_21_initial, c_22_initial, s_21_initial, s_22_initial, k_re, k_im, A,  B, C]';
delta_x    =  0.5 .* x;
zero_delta = zeros(10,1);

A = zeros(3,10);

T_g = coefficient_T_g .* [((1/sqrt(3)) .* c_20_initial - c_22_initial),        -s_22_initial,          -c_21_initial
                           -s_22_initial,              ((1/sqrt(3)) .* c_20_initial + c_22_initial),   -s_21_initial
                          -c_21_initial,                   -s_21_initial,               -(2/sqrt(3)) .* c_20_initial];
T_g = T_g + eye(3).*tr/3;        

t = 1;
stop = 20; % maximum number of iterations
threshold = 0.1; % 10% threshold
minimum   = 1e-3;
while t <= stop
    
    dh = h(:,t+1) - h(:,t);
    
    [T_g_plus_delta_x, w_plus_delta_x] = calculate_w(x, dt, w_in, T_g, h(:,t), r_sun, r_moon, dh, coefficient_F,coefficient_T_g, coefficient_T_r, GM_sun, GM_moon, delta_x);
    [T_g, w_0]                         = calculate_w(x, dt, w_in, T_g, h(:,t), r_sun, r_moon, dh, coefficient_F,coefficient_T_g, coefficient_T_r, GM_sun, GM_moon, zero_delta);
      
%     diff = (w_plus_delta_x(1)-w_0(1))/delta_x(1)
      
    A = getA(w_plus_delta_x, w_0, delta_x, minimum);
    
    l_0     = w_0;
    l       = reference(:,t);
    delta_l = l - l_0;
    
    delta_x_dach = ((A' * A) \ A') * delta_l;
    x_dach       = x + delta_x_dach;
    
    [T_g, w_new] = calculate_w(x_dach, dt, w_in, T_g, h(:,t), r_sun, r_moon, dh, coefficient_F,coefficient_T_g, coefficient_T_r, GM_sun, GM_moon, delta_x);
    
    % relativ difference to the previous iteration
    diff = abs((w_new-w_0)./w_0);
    % stop if difference to preious is below threshold value (default: 10%)
    if (max(diff) < threshold)
        break;
    end
      
    x    = x_dach;
    w_in = w_new;
    
    if isnan(w_in)
        disp(["Iteration: ", num2str(t)])
        disp("w_in is NAN: ")
        break
    end
    
    t    = t + 1;
end






