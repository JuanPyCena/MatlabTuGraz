clear all
close all
clc
format long
diary('logOutput.txt')
diary on
tic

max_iter = 3;
iter     = 1;
threshold_value = 1e-9;

G       = (6.674e-11) * 3600 * 3600;          % [m^3/(kg* h^2)]
GM_sun  = (1.32712442076e20) * 3600 * 3600;   % [m^3/ h^2]
GM_moon = (4.9027779e12) * 3600 * 3600 ;      % [m^3/ h^2]
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

%% Einlesen
AAM_FILE = 'ESMGFZ_AAM_v1.0_03h_2004.asc';
HAM_FILE = 'ESMGFZ_HAM_v1.2_24h_2004.asc';
OAM_FILE = 'ESMGFZ_OAM_v1.0_03h_2004.asc';

[timespan, step, h, grav_potent,r_moon,r_sun, reference] = read_data(AAM_FILE);
c_20 = grav_potent(1, :);
c_21 = grav_potent(2, :);
s_21 = grav_potent(3, :);
c_22 = grav_potent(4, :);
s_22 = grav_potent(5, :);

w_initial = reference(:,1) .*3600;

%% delta_X Vektor
%          w_initial_x,  w_initial_y,  w_initial_z,  k_re, k_im, tr  
delta_x = [1e-8,         1e-8,           1e-8,       1e-8, 1e-8, 1e-8];

while iter <= max_iter
    disp(['Iteration: ', num2str(iter)]);
    
    %% F(x)
    w_initial_x = [w_initial(1) + delta_x(1); w_initial(2);              w_initial(3)];
    w_initial_y = [w_initial(1);              w_initial(2) + delta_x(2); w_initial(3)];
    w_initial_z = [w_initial(1);              w_initial(2);              w_initial(3) + delta_x(3)];

    omega_delta_x = f_omega(w_initial_x, ...
                            r_sun, r_moon,...
                            c_20, c_21, c_22, s_21, s_22, h,...
                            coefficient_T_g, coefficient_T_r, coefficient_F,...
                            GM_sun, GM_moon, k_re, k_im, A, B, C, tr, timespan);

    omega_delta_y = f_omega(w_initial_y, ...
                            r_sun, r_moon,...
                            c_20, c_21, c_22, s_21, s_22, h,...
                            coefficient_T_g, coefficient_T_r, coefficient_F,...
                            GM_sun, GM_moon, k_re, k_im, A, B, C, tr, timespan);

    omega_delta_z = f_omega(w_initial_z , ...
                            r_sun, r_moon,...
                            c_20, c_21, c_22, s_21, s_22, h,...
                            coefficient_T_g, coefficient_T_r, coefficient_F,...
                            GM_sun, GM_moon, k_re, k_im, A, B, C, tr, timespan);
    
    omega_delta_k_re = f_omega(w_initial, ...
                               r_sun, r_moon,...
                               c_20, c_21, c_22, s_21, s_22, h,...
                               coefficient_T_g, coefficient_T_r, coefficient_F,...
                               GM_sun, GM_moon, k_re  + delta_x(4), k_im, A, B, C, tr, timespan);

    omega_delta_k_im = f_omega(w_initial, ...
                               r_sun, r_moon,...
                               c_20, c_21, c_22, s_21, s_22, h,...
                               coefficient_T_g, coefficient_T_r, coefficient_F,...
                               GM_sun, GM_moon, k_re, k_im  + delta_x(5), A, B, C, tr, timespan);
                         
    omega_delta_k_tr = f_omega(w_initial, ...
                               r_sun, r_moon,...
                               c_20, c_21, c_22, s_21, s_22, h,...
                               coefficient_T_g, coefficient_T_r, coefficient_F,...
                               GM_sun, GM_moon, k_re, k_im, A, B, C, tr  * (1+ delta_x(6)), timespan);
                         
    omega_0     = f_omega(w_initial,...
                          r_sun, r_moon,...
                          c_20, c_21, c_22, s_21, s_22, h,...
                          coefficient_T_g, coefficient_T_r, coefficient_F,...
                          GM_sun, GM_moon, k_re, k_im, A, B, C, tr, timespan);

    %% Prepare omega for matrix
    omega_delta_mat      = [];
    omega_delta_mat_x    = [];
    omega_delta_mat_y    = [];
    omega_delta_mat_z    = [];
    omega_delta_mat_k_re = [];
    omega_delta_mat_k_im = [];
    omega_delta_mat_tr   = [];
    omega_0_mat          = [];

    % Create big column vector of stacked omega_delta vectors
    for i = 1:3:numel(omega_delta_x)
        omega_delta_mat_x    = [omega_delta_mat_x, omega_delta_x(i), omega_delta_x(i+1), omega_delta_x(i+2)];
        omega_delta_mat_y    = [omega_delta_mat_y, omega_delta_y(i), omega_delta_y(i+1), omega_delta_y(i+2)];
        omega_delta_mat_z    = [omega_delta_mat_z, omega_delta_z(i), omega_delta_z(i+1), omega_delta_z(i+2)];
        omega_delta_mat_k_re = [omega_delta_mat_k_re, omega_delta_k_re(i), omega_delta_k_re(i+1), omega_delta_k_re(i+2)];
        omega_delta_mat_k_im = [omega_delta_mat_k_im, omega_delta_k_im(i), omega_delta_k_im(i+1), omega_delta_k_im(i+2)];
        omega_delta_mat_tr   = [omega_delta_mat_tr, omega_delta_k_tr(i), omega_delta_k_tr(i+1), omega_delta_k_tr(i+2)];
    end
    omega_delta_mat_x    = omega_delta_mat_x';
    omega_delta_mat_y    = omega_delta_mat_y';
    omega_delta_mat_z    = omega_delta_mat_z';
    omega_delta_mat_k_re = omega_delta_mat_k_re';
    omega_delta_mat_k_im = omega_delta_mat_k_im';
    omega_delta_mat_tr   = omega_delta_mat_tr';
    omega_delta_mat      = [omega_delta_mat_x, omega_delta_mat_y, ...
                            omega_delta_mat_z, omega_delta_mat_k_re, ...
                            omega_delta_mat_k_im, omega_delta_mat_tr];

    % Create big column vector of stacked omega_0 vectors
    for i = 1:3:numel(omega_0)
        omega_0_mat = [omega_0_mat, omega_0(i), omega_0(i+1), omega_0(i+2)];
    end
    omega_0_mat = omega_0_mat';

    %% A_mat
    A_mat = zeros(length(omega_delta_mat), numel(delta_x));

    % Amatrix(:,i) = (f(x0 + dx) - f(x0) / dx) 
    for i = 1:numel(delta_x)
       A_mat(:,i) = (omega_delta_mat(:,i) - omega_0_mat) ./ delta_x(i); 
    end

    %% Reduzierte Beobachtungen
    l = [];
    % Bring the reference Data into a big column vector like omega_0
    for i = 1:3:numel(omega_delta_x)
        l = [l, reference(i) .*3600,  reference(i+1) .*3600, reference(i+2) .*3600];
    end
    l = l';

    l_0 = omega_0_mat;

    delta_l = l - l_0;

    %% Schätzung der Loesung
    pseudo_inverse = (A_mat' * A_mat) \ A_mat';
    delta_x_dach   =  pseudo_inverse * delta_l;

    %% Korrigierte Lösung
    w_initial_corrected = [w_initial(1) + delta_x_dach(1);
                           w_initial(2) + delta_x_dach(2);
                           w_initial(3) + delta_x_dach(3)];
    omega_corrected = f_omega(w_initial_corrected,...
                              r_sun, r_moon,...
                              c_20, c_21, c_22, s_21, s_22, h,...
                              coefficient_T_g, coefficient_T_r, coefficient_F,...
                              GM_sun, GM_moon, k_re + delta_x_dach(4),...
                              k_im + delta_x_dach(5), A, B, C,...
                              tr + delta_x_dach(6), timespan);             
    
    %% Abbruchsbedingung 
    diff_corrected = omega_corrected - reference(:,1:length(omega_corrected)).*3600;
    diff_corrected_mat{iter} = diff_corrected;
    abbruch = abs(max(max(diff_corrected)));
    if (abbruch <= threshold_value)
        disp(['Difference smaller than: ', num2str(threshold_value)]);
        break;
    end
    
    %% Delta_X Vektor überschreiben, x_dach überschreiben und iter erhoehen.
    delta_x = delta_x_dach;
    iter    = iter + 1;
end

%% Plot
xp           = (R/omega_N) .* omega_corrected(1,:);
yp           = (R/omega_N) .* omega_corrected(2,:);

xp_estimated = (R/omega_N) .* omega_0(1,:);
yp_estimated = (R/omega_N) .* omega_0(2,:);

xp_reference = (R/omega_N) .* reference(1,1:length(omega_corrected)).*3600;
yp_reference = (R/omega_N) .* reference(2,1:length(omega_corrected)).*3600;

figure(1)
hold on
plot(xp,yp)
plot(xp_estimated, yp_estimated)
plot(xp_reference,yp_reference)
hold off
title('Polar Motion at Earth Surface')
legend('omega corrected','estimated omega', 'reference data')
ylabel('y[m]')
xlabel('x[m]')
axis equal
savefig('PolarPlot.fig')

toc
diary off













