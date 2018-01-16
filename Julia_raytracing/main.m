%% UE Methods of Modeling and Simulation
% Ray tracing through the atmosphere
% Sammer Julia, 01204903
% 5.1.2018

clear all 
close all
clc

%% Beispiel 1 - Polarplot Euler / Runge-Kutta
% Euler(phi0)
[e_theta0, r] = F_Euler(0);
[e_theta30, r] = F_Euler(30);
[e_theta60, r] = F_Euler(60);
[e_theta90, r] = F_Euler(89.9);

thetaEuler = [e_theta0; e_theta30; e_theta60; e_theta90];

% Runge-Kutta(phi0)
[rk_theta0] = F_RungeKutta(0);
[rk_theta30] = F_RungeKutta(30);
[rk_theta60] = F_RungeKutta(60);
[rk_theta90] = F_RungeKutta(89.9);

thetaRungeKutta = [rk_theta0; rk_theta30; rk_theta60; rk_theta90];

% Plot Euler
figure
polarplot(thetaEuler(:,2:end),r-6300)
ax = gca;
ax.ThetaLim = [0 90];

%Plot Runge-Kutta
figure
polarplot(thetaRungeKutta(:,2:end),r-6300)
ax = gca;
ax.ThetaLim = [0 90];

%Differenzenplot
figure
hold on
plot(thetaEuler(1,:) - thetaRungeKutta(1,:))
plot(thetaEuler(2,:) - thetaRungeKutta(2,:))
plot(thetaEuler(3,:) - thetaRungeKutta(3,:))
plot(thetaEuler(4,:) - thetaRungeKutta(4,:))












