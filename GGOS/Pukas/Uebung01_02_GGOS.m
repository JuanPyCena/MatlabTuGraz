%% clear memory:
clear all
close all
clc
%% Variable Declaration:

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

wx = 1.336419911142783967e-11;
wy = -5.587677166668280574e-11;
wz = 7.2921151467064e-5;
w = [wx;wy;wz].*3600;
%% Read Data:

AAM_FILE = 'ESMGFZ_AAM_v1.0_03h_2004.asc';
HAM_FILE = 'ESMGFZ_HAM_v1.2_24h_2004.asc';
OAM_FILE = 'ESMGFZ_OAM_v1.0_03h_2004.asc';

delimiterIn = ' ';
DATA_motionMomentumAAM = importdata(AAM_FILE, delimiterIn);
DATA_motionMomentumHAM = importdata(HAM_FILE, delimiterIn);
DATA_motionMomentumOAM = importdata(OAM_FILE, delimiterIn);
DATA_gravPotential  = importdata('potentialCoefficientsAOHIS.txt', delimiterIn);
DATA_moon           = importdata('moon.txt', delimiterIn);
DATA_sun            = importdata('sun.txt', delimiterIn);
DATA_reference      = importdata('earthRotationVector.txt', delimiterIn);

% allocating the timespan
timespan = length(DATA_gravPotential(:,1));
  
% defining the timestep by checking the fileName
if ~isempty(strfind(OAM_FILE, '_24h_'))
    step = 24
elseif ~isempty(strfind(OAM_FILE, '_03h_'))
    step = 3
else
    step = 1
end
    
% allocating the motionMomentum
h_readAAM(1,:) = DATA_motionMomentumAAM(:,9)  .* (omega_N * (C - A) / (1.610));
h_readAAM(2,:) = DATA_motionMomentumAAM(:,10) .* (omega_N * (C - A) / (1.610));
h_readAAM(3,:) = DATA_motionMomentumAAM(:,11) .* (omega_N *  C / (1.125));

h_readHAM(1,:) = DATA_motionMomentumHAM(:,9)  .* (omega_N * (C - A) / (1.610));
h_readHAM(2,:) = DATA_motionMomentumHAM(:,10) .* (omega_N * (C - A) / (1.610));
h_readHAM(3,:) = DATA_motionMomentumHAM(:,11) .* (omega_N *  C / (1.125));

h_readOAM(1,:) = DATA_motionMomentumOAM(:,9)  .* (omega_N * (C - A) / (1.610));
h_readOAM(2,:) = DATA_motionMomentumOAM(:,10) .* (omega_N * (C - A) / (1.610));
h_readOAM(3,:) = DATA_motionMomentumOAM(:,11) .* (omega_N *  C / (1.125));

hHAM(1,:) = (interp1(h_readHAM(1,:),1:1/8:366.8750))';
hHAM(2,:) = (interp1(h_readHAM(2,:),1:1/8:366.8750))';
hHAM(3,:) = (interp1(h_readHAM(3,:),1:1/8:366.8750))';

h_read = h_readAAM  + h_readOAM + hHAM;

h(1,:) = (interp1(h_read(1,:),1:1/step:length(h_read)))';
h(2,:) = (interp1(h_read(2,:),1:1/step:length(h_read)))';
h(3,:) = (interp1(h_read(3,:),1:1/step:length(h_read)))';

% allocating the gravitational potential
c_20 = DATA_gravPotential(1:timespan,2);
c_21 = DATA_gravPotential(1:timespan,3);
s_21 = DATA_gravPotential(1:timespan,4);
c_22 = DATA_gravPotential(1:timespan,5);
s_22 = DATA_gravPotential(1:timespan,6);
 
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

%% pre-caclulations (Momentum):

dist_sun   = sqrt(r_sun(1,:).^2 + r_sun(2,:).^2 + r_sun(3,:).^2);
dist_moon   = sqrt(r_moon(1,:).^2 + r_moon(2,:).^2 + r_moon(3,:).^2);

M_sun(1,:) = (((C-B) .* r_sun(2,:) .* r_sun(3,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);
M_sun(2,:) = (((A-C) .* r_sun(1,:) .* r_sun(3,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);
M_sun(3,:) = (((B-A) .* r_sun(1,:) .* r_sun(2,:)) ./ (dist_sun.^5)) .* (3 * GM_sun);

M_moon(1,:) = (((C-B) .* r_moon(2,:) .* r_moon(3,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
M_moon(2,:) = (((A-C) .* r_moon(1,:) .* r_moon(3,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);
M_moon(3,:) = (((B-A) .* r_moon(1,:) .* r_moon(2,:)) ./ (dist_moon.^5)) .* (3 * GM_moon);

M = M_sun + M_moon;

%% pre-calculations  (differential variables - TG, d_TG, dh):
I = eye(3);

for i = 1:timespan
    T_g{i}  =  tr/3.*I + (sqrt(5/3)*Mass*R^2) .* [1/sqrt(3)*c_20(i)-c_22(i), -s_22(i), -c_21(i);
                                                    -s_22(i), 1/sqrt(3)*c_20(i)+c_22(i), -s_21(i);
                                                    -c_21(i), -s_21(i), -2/sqrt(3)*c_20(i)];
                                                                                         
end
for time = 1:timespan-1
    dT_g{time} = T_g{time + 1} - T_g{time};
end

for time = 1:(length(h)-1)
    dh{time} = h(:,time + 1) - h(:,time);
end

%% Differential Equation:

for t = 1:24*365
    
    T_r = coefficient_T_r .* [0,     0,      k_re * w(1) + k_im * w(2);
                              0,     0,      k_re * w(2) - k_im * w(1);
            k_re * w(1) + k_im * w(2),     k_re * w(2) - k_im * w(1),         0];
        

    F = T_g{t} + T_r + coefficient_F .* [k_re, k_im, 0;
                                        -k_im, k_re, 0;
                                         0   ,  0  , 0];
    
    w = w + (F\(M(:,t) - (dT_g{t} * w) - cross(w, (T_g{t} + T_r) * w) - cross(w,h(:,t)) - dh{t})); 
    
    result{t} = w;
    
    disp ('Tag: ')
    disp (t/24) 
end

%% Uebung 2 Least Square:

for t = 1:24*365
    
    T_r = coefficient_T_r .* [0,     0,      k_re * w(1) + k_im * w(2);
                              0,     0,      k_re * w(2) - k_im * w(1);
            k_re * w(1) + k_im * w(2),     k_re * w(2) - k_im * w(1),         0];
        

    F = T_g{t} + T_r + coefficient_F .* [k_re, k_im, 0;
                                        -k_im, k_re, 0;
                                         0   ,  0  , 0];
    
    w = w + (F\(M(:,t) - (dT_g{t} * w) - cross(w, (T_g{t} + T_r) * w) - cross(w,h(:,t)) - dh{t})); 
    
    xp1        = (R/omega_N) .* w(1);
    yp1        = (R/omega_N) .* w(2);

    xp_reference = (R/omega_N) .* reference(1,t) .* 3600;
    yp_reference = (R/omega_N) .* reference(2,t) .* 3600;
    
    RSS1_x = ((reference(1,t)-w(1))*(xp_reference - xp1))/((reference(1,t)-w(1))^2);
    RSS1_y = ((reference(2,t)-w(2))*(yp_reference - yp1))/((reference(2,t)-w(2))^2);
    
    RSS0_x = xp1 - RSS1_x*w(1);
    RSS0_y = yp1 - RSS1_y*w(2);
    
    xp{t} = xp1 - RSS0_x;
    yp{t} = yp1 - RSS0_y;
    
    result{t} = w;
    
    disp ('Tag: ')
    disp (t/24) 
end


%% Plotting: 

result_m = cell2mat(result);
% 
xp        = (R/omega_N) .* result_m(1,:);
yp        = (R/omega_N) .* result_m(2,:);

% xp_reference = (R/omega_N) .* reference(1,:) .* 3600;
% yp_reference = (R/omega_N) .* reference(2,:) .* 3600;

xp_m = cell2mat(xp);
yp_m = cell2mat(yp);

% delta_LOD = 86400/omega_N .* (omega_N - result_m(3,:)./3600);
% LOD = 8.6376e+04 - delta_LOD;
% LOD = LOD .* 1000;
figure(1)
plot(xp,yp)
hold on
plot(xp_reference(1:8760),yp_reference(1:8760))
title('Polar Motion at Earth Surface')
ylabel('y[m]')
xlabel('x[m]')
saveas(gcf,'PolarMotion.png')
% hold off
% 
% figure(2)
% plot(delta_LOD)
% title('Length of day')
% ylabel('Length of Day [s]')
% xlabel('Time [h]')
% 
% figure(3)
% subplot(3,1,1)
% plot(result_m(1,:)./3600)
% title('Angular velocity - w_x')
% ylabel('w_x [rad/s]')
% xlabel('Time [h]')
% subplot(3,1,2)
% plot(result_m(2,:)./3600)
% title('Angular velocity - w_y')
% ylabel('w_y [rad/s]')
% xlabel('Time [h]')
% subplot(3,1,3)
% plot(result_m(3,:)./3600)
% title('Angular velocity - w_z')
% ylabel('w_z [rad/s]')
% xlabel('Time [h]')
% saveas(gcf,'AngularVelocity.png')
% 
% figure(4)
% plot(LOD)
% title('Delta Length of day')
% ylabel('Length of Day [ms]')
% xlabel('Time [h]')
