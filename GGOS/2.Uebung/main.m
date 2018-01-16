clear all;
close all;
clc
%% General Values

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


%% Calculation
dt = 1;
counter = 1;
w_in = w_initial.*3600;


x          = [c_20_initial, c_21_initial, c_22_initial, s_21_initial, s_22_initial, k_re, k_im, A,  B, C]';
delta_x    = [1,1,1,1,1,1,1,1,1,1]';
zero_delta = zeros(10,1);

A = zeros(3,10);

T_g = coefficient_T_g .* [((1/sqrt(3)) .* c_20_initial - c_22_initial) + tr/3,        -s_22_initial,          -c_21_initial
                           -s_22_initial,              ((1/sqrt(3)) .* c_20_initial + c_22_initial) + tr/3,   -s_21_initial
                          -c_21_initial,                   -s_21_initial,               -(2/sqrt(3)) .* c_20_initial + tr/3];
                      
for t = 1:dt:1
    
    dh = h(:,t+1) - h(:,t);
    
   [T_g_plus_delta_x, w_plus_delta_x] = calculate_w(x, dt, w_in, T_g, h(:,t), r_sun, r_moon, dh, coefficient_F,coefficient_T_g, coefficient_T_r, GM_sun, GM_moon, delta_x);
   [T_g, w_0]                         = calculate_w(x, dt, w_in, T_g, h(:,t), r_sun, r_moon, dh, coefficient_F,coefficient_T_g, coefficient_T_r, GM_sun, GM_moon, zero_delta);
   
   %Elements for A
   %w_x derivatives for parameters
   A(1,1)  = (w_plus_delta_x(1) - w_0(1))/delta_x(1); 
   A(1,2)  = (w_plus_delta_x(1) - w_0(1))/delta_x(2); 
   A(1,3)  = (w_plus_delta_x(1) - w_0(1))/delta_x(3); 
   A(1,4)  = (w_plus_delta_x(1) - w_0(1))/delta_x(4); 
   A(1,5)  = (w_plus_delta_x(1) - w_0(1))/delta_x(5); 
   A(1,6)  = (w_plus_delta_x(1) - w_0(1))/delta_x(6); 
   A(1,7)  = (w_plus_delta_x(1) - w_0(1))/delta_x(7); 
   A(1,8)  = (w_plus_delta_x(1) - w_0(1))/delta_x(8);
   A(1,9)  = (w_plus_delta_x(1) - w_0(1))/delta_x(9); 
   A(1,10) = (w_plus_delta_x(1) - w_0(1))/delta_x(10); 
   
   %w_y derivatives for parameters
   A(2,1)  = (w_plus_delta_x(2) - w_0(2))/delta_x(1); 
   A(2,2)  = (w_plus_delta_x(2) - w_0(2))/delta_x(2); 
   A(2,3)  = (w_plus_delta_x(2) - w_0(2))/delta_x(3); 
   A(2,4)  = (w_plus_delta_x(2) - w_0(2))/delta_x(4); 
   A(2,5)  = (w_plus_delta_x(2) - w_0(2))/delta_x(5); 
   A(2,6)  = (w_plus_delta_x(2) - w_0(2))/delta_x(6); 
   A(2,7)  = (w_plus_delta_x(2) - w_0(2))/delta_x(7); 
   A(2,8)  = (w_plus_delta_x(2) - w_0(2))/delta_x(8);
   A(2,9)  = (w_plus_delta_x(2) - w_0(2))/delta_x(9); 
   A(2,10) = (w_plus_delta_x(2) - w_0(2))/delta_x(10); 
   
   %w_z derivatives for parameters
   A(3,1)  = (w_plus_delta_x(3) - w_0(3))/delta_x(1); 
   A(3,2)  = (w_plus_delta_x(3) - w_0(3))/delta_x(2); 
   A(3,3)  = (w_plus_delta_x(3) - w_0(3))/delta_x(3); 
   A(3,4)  = (w_plus_delta_x(3) - w_0(3))/delta_x(4); 
   A(3,5)  = (w_plus_delta_x(3) - w_0(3))/delta_x(5); 
   A(3,6)  = (w_plus_delta_x(3) - w_0(3))/delta_x(6); 
   A(3,7)  = (w_plus_delta_x(3) - w_0(3))/delta_x(7); 
   A(3,8)  = (w_plus_delta_x(3) - w_0(3))/delta_x(8);
   A(3,9)  = (w_plus_delta_x(3) - w_0(3))/delta_x(9); 
   A(3,10) = (w_plus_delta_x(3) - w_0(3))/delta_x(10); 
   
   disp(A * delta_x)
end






