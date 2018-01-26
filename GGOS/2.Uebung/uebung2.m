clear all
close all
clc
format long
diary('logOutput.txt')
diary on
tic

year_begin         = 4;
year_end           = 15;
numOfDays          = 365 * (year_end - year_begin);
max_iter           = 420;
threshold_value    = 1e-8;
threshold_relative = 1e-12;

variations_in_calculation = ["Only_omega0_is_altered"
                             "Omega0_and_the_love_numbers_are_altered" 
                             "Omega0_and_the_love_numbers_and_the_trace_are_altered"];
                         
fileNames = ["onlyOmega.txt"
             "omegaAndLove.txt"
             "omegaLoveAndtrace.txt"];

%% General Values
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

%% Calibrating and calculating all omegas
for variation = 1:numel(variations_in_calculation)
    abbruch_m          = [];
    abbruch_rel_m      = [];
    iter               = 1;

    %% Einlesen
    AAM_FILE = 'ESMGFZ_AAM_v1.0_03h_20';

    [timespan, step, h, grav_potent,r_moon,r_sun, reference] = read_data(AAM_FILE, year_begin, year_end);
    c_20 = grav_potent(1, :);
    c_21 = grav_potent(2, :);
    s_21 = grav_potent(3, :);
    c_22 = grav_potent(4, :);
    s_22 = grav_potent(5, :);

    w_initial = reference(:,1) .*3600;

    %% delta_X Vektor
    %          w_initial_x,  w_initial_y,  w_initial_z,  k_re, k_im, tr  
    delta_val = 1;
    
    % Only omega0 is altered
    if variation == 1
        delta_x   = [delta_val,       delta_val,    delta_val];
        x_vec     = [w_initial(1), w_initial(2), w_initial(3)]'; 
    % Omega0 and the love numbers are altered    
    elseif variation == 2
        delta_x   = [delta_val,       delta_val,    delta_val, delta_val, delta_val];
        x_vec     = [w_initial(1), w_initial(2), w_initial(3), k_re     , k_im]';
    % Omega0, the love numbers and the trace are altered    
    elseif variation == 3        
        delta_x   = [delta_val,       delta_val,    delta_val, delta_val, delta_val, 1e38];
        x_vec     = [w_initial(1), w_initial(2), w_initial(3), k_re     , k_im     , tr]';
    end

    %% Iterieren
    while iter <= max_iter
        disp(['Iteration: ', num2str(iter)]);

        %% X_vektor assignen
        w_initial = [x_vec(1); x_vec(2); x_vec(3)];
        
        % Omega0 and the love numbers are altered    
        if variation >= 2
            k_re      = [x_vec(4)];
            k_im      = [x_vec(5)];
        end
        
        % Omega0, the love numbers and the trace are altered   
        if variation == 3
            tr        = [x_vec(6)];
        end

        %% F(x)
        w_initial_x = [w_initial(1) + delta_x(1); w_initial(2);              w_initial(3)];
        w_initial_y = [w_initial(1);              w_initial(2) + delta_x(2); w_initial(3)];
        w_initial_z = [w_initial(1);              w_initial(2);              w_initial(3) + delta_x(3)];

        omega_delta_x = f_omega(w_initial_x, ...
                                r_sun, r_moon,...
                                c_20, c_21, c_22, s_21, s_22, h,...
                                coefficient_T_g, coefficient_T_r, coefficient_F,...
                                GM_sun, GM_moon, k_re, k_im, A, B, C, tr, timespan, numOfDays);

        omega_delta_y = f_omega(w_initial_y, ...
                                r_sun, r_moon,...
                                c_20, c_21, c_22, s_21, s_22, h,...
                                coefficient_T_g, coefficient_T_r, coefficient_F,...
                                GM_sun, GM_moon, k_re, k_im, A, B, C, tr, timespan, numOfDays);
        
        omega_delta_z = f_omega(w_initial_z , ...
                                r_sun, r_moon,...
                                c_20, c_21, c_22, s_21, s_22, h,...
                                coefficient_T_g, coefficient_T_r, coefficient_F,...
                                GM_sun, GM_moon, k_re, k_im, A, B, C, tr, timespan, numOfDays);

        % Omega0 and the love numbers are altered    
        if variation >= 2
            omega_delta_k_re = f_omega(w_initial, ...
                                       r_sun, r_moon,...
                                       c_20, c_21, c_22, s_21, s_22, h,...
                                       coefficient_T_g, coefficient_T_r, coefficient_F,...
                                       GM_sun, GM_moon, k_re  + delta_x(4), k_im, A, B, C, tr, timespan, numOfDays);

            omega_delta_k_im = f_omega(w_initial, ...
                                       r_sun, r_moon,...
                                       c_20, c_21, c_22, s_21, s_22, h,...
                                       coefficient_T_g, coefficient_T_r, coefficient_F,...
                                       GM_sun, GM_moon, k_re, k_im  + delta_x(5), A, B, C, tr, timespan, numOfDays);
        end
        
        % Omega0, the love numbers and the trace are altered
        if variation == 3
            omega_delta_k_tr = f_omega(w_initial, ...
                                       r_sun, r_moon,...
                                       c_20, c_21, c_22, s_21, s_22, h,...
                                       coefficient_T_g, coefficient_T_r, coefficient_F,...
                                       GM_sun, GM_moon, k_re, k_im, A, B, C, tr + delta_x(6), timespan, numOfDays);
        end
        
        omega_0     = f_omega(w_initial,...
                              r_sun, r_moon,...
                              c_20, c_21, c_22, s_21, s_22, h,...
                              coefficient_T_g, coefficient_T_r, coefficient_F,...
                              GM_sun, GM_moon, k_re, k_im, A, B, C, tr, timespan, numOfDays);

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
            
            % Omega0 and the love numbers are altered    
            if variation >= 2
                omega_delta_mat_k_re = [omega_delta_mat_k_re, omega_delta_k_re(i), omega_delta_k_re(i+1), omega_delta_k_re(i+2)];
                omega_delta_mat_k_im = [omega_delta_mat_k_im, omega_delta_k_im(i), omega_delta_k_im(i+1), omega_delta_k_im(i+2)];
            end
            
            % Omega0, the love numbers and the trace are altered
            if variation == 3
                omega_delta_mat_tr   = [omega_delta_mat_tr, omega_delta_k_tr(i), omega_delta_k_tr(i+1), omega_delta_k_tr(i+2)];
            end
        end

        omega_delta_mat_x    = omega_delta_mat_x';
        omega_delta_mat_y    = omega_delta_mat_y';
        omega_delta_mat_z    = omega_delta_mat_z';
        omega_delta_mat_k_re = omega_delta_mat_k_re';
        omega_delta_mat_k_im = omega_delta_mat_k_im';
        omega_delta_mat_tr   = omega_delta_mat_tr';
        
        % Only omega0 is altered
        if variation == 1
            omega_delta_mat  = [omega_delta_mat_x, ...
                                omega_delta_mat_y, ...
                                omega_delta_mat_z];
                            
        % Omega0 and the love numbers are altered    
        elseif variation == 2
            omega_delta_mat  = [omega_delta_mat_x, ...
                                omega_delta_mat_y, ...
                                omega_delta_mat_z, ...
                                omega_delta_mat_k_re, ...
                                omega_delta_mat_k_im];
                            
        % Omega0, the love numbers and the trace are altered    
        elseif variation == 3        
            omega_delta_mat  = [omega_delta_mat_x, ...
                                omega_delta_mat_y, ...
                                omega_delta_mat_z, ...
                                omega_delta_mat_k_re, ...
                                omega_delta_mat_k_im, ...
                                omega_delta_mat_tr];
        end

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
        
        if any(isnan(delta_x_dach))
            disp('An Element of delta_x_dach is NAN! Exiting before anything bad happens');
            break;
        end

        %% Korrigierte Lösung
        w_initial_corrected = [w_initial(1) + delta_x_dach(1);
                               w_initial(2) + delta_x_dach(2);
                               w_initial(3) + delta_x_dach(3)];
                           
        % Only omega0 is altered            
        if variation == 1
                omega_corrected = f_omega(w_initial_corrected,...
                                  r_sun, r_moon,...
                                  c_20, c_21, c_22, s_21, s_22, h,...
                                  coefficient_T_g, coefficient_T_r, coefficient_F,...
                                  GM_sun, GM_moon, ...
                                  k_re,...
                                  k_im,...
                                  A, B, C,...
                                  tr,...
                                  timespan, numOfDays);  
        % Omega0 and the love numbers are altered    
        elseif variation == 2
                omega_corrected = f_omega(w_initial_corrected,...
                                  r_sun, r_moon,...
                                  c_20, c_21, c_22, s_21, s_22, h,...
                                  coefficient_T_g, coefficient_T_r, coefficient_F,...
                                  GM_sun, GM_moon, ...
                                  k_re + delta_x_dach(4),...
                                  k_im + delta_x_dach(5),...
                                  A, B, C,...
                                  tr ,...
                                  timespan, numOfDays);  
        % Omega0, the love numbers and the trace are altered    
        elseif variation == 3        
                omega_corrected = f_omega(w_initial_corrected,...
                                  r_sun, r_moon,...
                                  c_20, c_21, c_22, s_21, s_22, h,...
                                  coefficient_T_g, coefficient_T_r, coefficient_F,...
                                  GM_sun, GM_moon, ...
                                  k_re + delta_x_dach(4),...
                                  k_im + delta_x_dach(5),...
                                  A, B, C,...
                                  tr  + delta_x(6),...
                                  timespan, numOfDays);    
        end
        
        x_dach = x_vec + delta_x_dach;                

        %% Abbruchsbedingung 
        diff_corrected           = omega_corrected - reference(:,1:length(omega_corrected)).*3600;
        diff_corrected_mat{iter} = diff_corrected;
        abbruch                  = abs(max(max(diff_corrected)));
        abbruch_m                = [abbruch_m, abbruch];

        if ((abbruch <= threshold_value))        
            disp(['Difference smaller than: ', num2str(threshold_relative)]);
            break;
        end

        if iter > 1
            abbruch_rel   = abs(abbruch_m(iter) - abbruch_m(iter-1));
            abbruch_rel_m = [abbruch_rel_m, abbruch_rel];
            if ((abbruch_rel <= threshold_relative))
                disp(['Difference between difference smaller than: ', num2str(threshold_relative)]);
                break;
            end
        end

        %% Delta_X Vektor überschreiben, x_dach überschreiben und iter erhoehen.
        delta_x = delta_x_dach;
        x_vec   = x_dach;
        iter    = iter + 1;
    end

    %% Berechnung mit kalibrierten Daten

    numOfDays = 365 * (year_end - year_begin +1);

    % Neu einlesen für letzte berechnung
    [timespan, step, h, grav_potent,r_moon,r_sun, reference] = read_data(AAM_FILE, year_begin, year_end);
    c_20 = grav_potent(1, :);
    c_21 = grav_potent(2, :);
    s_21 = grav_potent(3, :);
    c_22 = grav_potent(4, :);
    s_22 = grav_potent(5, :);

    w_initial = [x_vec(1); x_vec(2); x_vec(3)];
        
   % Omega0 and the love numbers are altered    
    if variation >= 2
        k_re      = [x_vec(4)];
        k_im      = [x_vec(5)];
    end

    % Omega0, the love numbers and the trace are altered   
    if variation == 3
        tr        = [x_vec(6)];
    end

    omega_calibrated    = f_omega(w_initial,...
                                  r_sun, r_moon,...
                                  c_20, c_21, c_22, s_21, s_22, h,...
                                  coefficient_T_g, coefficient_T_r, coefficient_F,...
                                  GM_sun, GM_moon, k_re, k_im, A, B, C, tr, timespan, numOfDays);
                              
     %% Kalibrierte Daten als .txt speichern
  
     T = table(omega_calibrated(1,:)', omega_calibrated(2,:)', omega_calibrated(3,:)');
     T.Properties.Description   = char(variations_in_calculation(variation));
     T.Properties.VariableNames = {'omega_x' 'omega_y' 'omega_z'};
     T.Properties.VariableUnits = {'rad/h' 'rad/h' 'rad/h'};
     writetable(T, char(fileNames(variation)))
end
toc
diary off

%% Read in Data as rad/h for Plot
omega_N = (7.2921151467064e-5) .* 3600;        % [rad/h]
R       = 6378136.6;                           % [m]

result_omega   = table2array(readtable(char(fileNames(1))))';
result_love    = table2array(readtable(char(fileNames(2))))';
result_trace   = table2array(readtable(char(fileNames(3))))';

DATA_reference      = importdata('earthRotationVector.txt', ' ');
reference_data(1,:) = DATA_reference(1:length(result_omega),2)' .* 3600;
reference_data(2,:) = DATA_reference(1:length(result_omega),3)' .* 3600;
reference_data(3,:) = DATA_reference(1:length(result_omega),4)' .* 3600;

%% Preparing all relevant data for plotting

xp_omega     = (R/omega_N) .* result_omega(1,:);
yp_omega     = (R/omega_N) .* result_omega(2,:);
xp_love      = (R/omega_N) .* result_love(1,:);
yp_love      = (R/omega_N) .* result_love(2,:);
xp_trace     = (R/omega_N) .* result_trace(1,:);
yp_trace     = (R/omega_N) .* result_trace(2,:);

xp_reference = (R/omega_N) .* reference_data(1,:);
yp_reference = (R/omega_N) .* reference_data(2,:);

delta_LOD_omega = 86400/omega_N .* (omega_N - result_omega(3,:)./3600);
delta_LOD_love  = 86400/omega_N .* (omega_N - result_love(3,:)./3600);
delta_LOD_trace = 86400/omega_N .* (omega_N - result_trace(3,:)./3600);
delta_LOD_ref   = 86400/omega_N .* (omega_N - reference_data(3,:)./3600);

LOD_omega = 8.6376e+04 - delta_LOD_omega;
LOD_love  = 8.6376e+04 - delta_LOD_love;
LOD_trace = 8.6376e+04 - delta_LOD_trace;
LOD_ref   = 8.6376e+04 - delta_LOD_ref;
LOD_omega = LOD_omega .* 1000;
LOD_love  = LOD_love .* 1000;
LOD_trace = LOD_trace .* 1000;
LOD_ref   = LOD_ref .* 1000;

%% Plotting
figure(1)
hold on
plot(xp_omega,yp_omega, 'LineWidth',0.1)
plot(xp_love,yp_love, 'LineWidth',0.1)
plot(xp_trace,yp_trace, 'LineWidth',0.1)
plot(xp_reference,yp_reference, 'LineWidth',0.1)
hold off
title('Polar Motion at Earth Surface')
legend('Only \omega_0', '\omega_0 and Love Numbers', '\omega_0, Love Numbers and trace' ,'reference data')
ylabel('y[m]')
xlabel('x[m]')
axis equal
savefig('CalibratedPolarPlot.fig')

figure(2)
subplot(3,1,1)
hold on
plot(result_omega(1,:)./3600)
plot(result_love(1,:)./3600)
plot(result_trace(1,:)./3600)
hold off
title('Angular velocity - w_x')
ylabel('w_x [rad/s]')
xlabel('Time [h]')
subplot(3,1,2)
hold on
plot(result_omega(2,:)./3600)
plot(result_love(2,:)./3600)
plot(result_trace(2,:)./3600)
hold off
title('Angular velocity - w_y')
ylabel('w_y [rad/s]')
xlabel('Time [h]')
subplot(3,1,3)
hold on
plot(result_omega(3,:)./3600)
plot(result_love(3,:)./3600)
plot(result_trace(3,:)./3600)
hold off
title('Angular velocity - w_z')
ylabel('w_z [rad/s]')
xlabel('Time [h]')

figure(3)
hold on
plot(LOD_omega)
plot(LOD_love)
plot(LOD_trace)
plot(LOD_ref)
hold off
legend('Only \omega_0', '\omega_0 and Love Numbers', '\omega_0, Love Numbers and trace' ,'reference data')
title('Delta Length of day')
ylabel('Length of Day [ms]')
xlabel('Time [h]')



